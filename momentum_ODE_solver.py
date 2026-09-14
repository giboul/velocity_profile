r"""1D momentum balance for flow over and through a porous bed, in FEniCSx.

Solves, for the velocity u(z) on z in [z0, z1]::

    0 = eps*rho*g*i + eps*fbar + (tau_v + tau_t + tau_d)' - rho*nu*eps'*u'

    tau_v = rho*nu*(eps*u)'                                       viscous
    tau_t = rho*eps*(kappa*Z*D)^2 * |u'|*u'                    turbulent
    tau_d = rho*(lam*d_p/(1-eps_b))*(1-eps)*eps*u*u'           dispersive
    D     = 1 - exp(-sqrt(Z*u/(nu*A_star^2)))                  van Driest damping

    eps*fbar = -A_E*(1-eps)^2/(eps*d_p^2)*rho*nu*u                 Ergun, linear
               -B_E*(1-eps)*eps/d_p*rho*|u|*u                  Ergun, quadratic

z points UP: z < 0 is inside the bed (eps -> eps_b), z > 0 is the free water
(eps -> 1).  Z(z), the "integral depth", is the height above the interface.
The bed is held at the deep-bed Ergun equilibrium u_SSL; the free surface
carries zero total stress, which is the natural BC of the weak form and so
needs no code at all.

|u'| and |u| appear instead of u'^2 and u^2 so that the turbulent and Ergun
terms stay dissipative for any Newton iterate, not only for those with
u, u' > 0.  On the physical solution they are identical to the squares.

Every function takes the parameters it uses as arguments.  Only the physical
constants -- the ones a run never varies -- carry a default, taken from the
UPPERCASE block below; the variables of the case and the solver settings are
required, so no call can quietly inherit the wrong bed, fluid or mesh.  The
block is still the single place those numbers live: `main` reads it and passes
them down.

The profiles also take `engine`, the module their maths goes through.  It
defaults to `numpy`, so calling them plainly evaluates them on arrays; the weak
form is the one place that passes `ufl` instead.

Run:  python fenics/fenics_attempt.py
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
import ufl
from dolfinx import fem, mesh
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI


plt.style.use("communications/hyporheic_flow/lhe.mplstyle")
# Consistent naming between the two engines, so `engine.<f>` works for both.
ufl.maximum = ufl.max_value
ufl.log = ufl.ln


# ===========================================================================
# Parameters -- the reference case, used as the functions' default arguments
# ===========================================================================
# Constants
G = 9.81               # gravity                              [m/s^2]
KAPPA = 0.41           # von Karman constant                  [-]
LAM = 0.015            # dispersive-mixing coefficient        [-]
A_STAR = 26.0          # van Driest constant                  [-]
A_E = 180.0            # Ergun linear (Kozeny-Carman)         [-]
B_E = 1.75             # Ergun quadratic (Forchheimer)        [-]

# Variables
RHO = 950.0            # density                              [kg/m^3]
NU = 3.1e-6            # kinematic viscosity                  [m^2/s]
SLOPE = 1.0e-2         # hydraulic gradient                   [-]
D_P = 8.0e-3           # grain diameter                       [m]
EPS_B = 0.38           # bulk (deep-bed) porosity             [-]
Z1 = 12e-3          # free surface                         [m]
Z0 = - Z1          # bottom of the bed                    [m]

# Parameters for solver
N_CELLS = 10_000       # 1D, so we can afford to be generous
GRADING = 1.9          # mesh clustering towards z = 0 (0 = uniform)
REG = 1.0e-12          # smooths |.| and sqrt(.) at the origin so that their
                       # derivatives stay finite for Newton


# # %%
# RHO = 1000.0            # density                              [kg/m^3]
# NU = 1.0e-6            # kinematic viscosity                  [m^2/s]
# SLOPE = 1.0e-2         # hydraulic gradient                   [-]
# D_P = 8.0e-3           # grain diameter                       [m]
# EPS_B = 0.2           # bulk (deep-bed) porosity             [-]
# Z1 = 15 * D_P          # free surface                         [m]
# Z0 = - Z1
# # %%

# ===========================================================================
# Derived quantities
# ===========================================================================
def interface_sharpness(d_p):
    """a [1/m]: eps goes eps_b -> 99.9% of 1 over roughly one grain diameter."""
    return -2.0 / d_p * np.log(1.0 / 0.994 - 1.0)


# ===========================================================================
# Vertical profiles
# ===========================================================================
def porosity(z, d_p, eps_b, engine=np):
    """eps(z), a smooth step from eps_b in the bed to 1 in the free water."""
    a = interface_sharpness(d_p)
    return eps_b + (1.0 - eps_b) * 0.5 * (1.0 + engine.tanh(0.5 * a * z))


def porosity_rate(z, d_p, eps_b, engine=np):
    """eps'(z), analytically: (1-eps_b)*a/4 * sech^2(a*z/2)."""
    a = interface_sharpness(d_p)
    return (1.0 - eps_b) * a / 4.0 * (1.0 - engine.tanh(0.5 * a * z)**2)

def integral_depth(z, d_p, engine=np):
    """Z(z) = softplus(a*z)/a: height above the interface, smoothed over 1/a.

    Written as max(x, 0) + log(1 + exp(-|x|)) rather than log(1 + exp(x)):
    the two are the same function, but the second overflows float64 at
    a*z = 709, i.e. z = 51.3*d_p.  The inf then reaches tau_t as inf - inf
    and SNES aborts on iteration 0 with a NaN residual (reason -4).
    """
    a = interface_sharpness(d_p)
    return engine.maximum(z, 0.0) + engine.log(1.0 + engine.exp(-abs(a * z))) / a


def mixing_length(z, u, d_p, nu, reg, kappa=KAPPA, A_star=A_STAR, engine=np):
    """Prandtl length kappa*Z with van Driest near-wall damping."""
    Z = integral_depth(z, d_p=d_p, engine=engine)
    u_pos = engine.maximum(u, 0.0)
    damping = 1.0 - engine.exp(-engine.sqrt(Z * u_pos / (nu * A_star**2) + reg))
    return kappa * (
        Z  # engine.maximum(z, 0)
        # + 2.5*d_p/30  # TODO : remove
    ) * damping


# ===========================================================================
# Stresses and body forces
# ===========================================================================
def stresses(z, u, du, d_p, eps_b, rho, nu, reg, kappa=KAPPA, lam=LAM,
             A_star=A_STAR, engine=np):
    """(tau_v, tau_t, tau_d) -- viscous, turbulent and dispersive."""
    eps = porosity(z, d_p=d_p, eps_b=eps_b, engine=engine)
    deps = porosity_rate(z, d_p=d_p, eps_b=eps_b, engine=engine)
    ell = mixing_length(z, u, d_p=d_p, nu=nu, reg=reg, kappa=kappa,
                        A_star=A_star, engine=engine)
    tau_v = rho * nu * (deps * u + eps * du)
    tau_t = rho * eps * ell**2 * engine.sqrt(du * du + reg) * du
    tau_d = rho * lam * d_p / (1.0 - eps_b) * (1.0 - eps) * eps * u * du
    return tau_v, tau_t, tau_d


def drag(e, u, d_p, rho, nu, reg, A_E=A_E, B_E=B_E, engine=np):
    """Ergun resistance eps*fbar: Kozeny-Carman plus Forchheimer."""
    return (
        - A_E * (1.0 - e)**2 / (e**2 * d_p**2) * rho * nu * u
        - B_E * (1.0 - e) / d_p * rho * engine.sqrt(u * u + reg) * u
    )


def source(z, u, du, slope, d_p, eps_b, rho, nu, reg, g=G, A_E=A_E, B_E=B_E,
           engine=np):
    """Every term of the strong form except tau', i.e. 0 = source(u) + tau'."""
    e = porosity(z, d_p=d_p, eps_b=eps_b, engine=engine)
    deps = porosity_rate(z, d_p=d_p, eps_b=eps_b, engine=engine)
    f = drag(e, u, d_p=d_p, rho=rho, nu=nu, reg=reg, A_E=A_E, B_E=B_E,
             engine=engine)
    return e * rho * g * slope + e * f - rho * nu * deps * du


def u_subsurface(slope, d_p, eps_b, rho, nu, g=G, A_E=A_E, B_E=B_E):
    """Deep-bed equilibrium: the positive root of eps*rho*g*i + eps*fbar = 0."""
    qa = - rho * B_E * eps_b * (1.0 - eps_b) / d_p
    qb = - rho * nu * A_E * (1.0 - eps_b)**2 / (eps_b * d_p**2)
    qc = + eps_b * rho * g * slope
    return (-qb - np.sqrt(qb**2 - 4.0 * qa * qc)) / (2.0 * qa)


# ===========================================================================
# Mesh
# ===========================================================================
def graded_mesh(z0, z1, n_cells, grading):
    """Interval [z0, z1], graded towards z = 0 where the porosity layer is.

    The map is monotone in zeta -- so the mesh stays valid whatever the vertex
    ordering -- and increasingly dense towards zeta = 0, which maps to z = 0.

    Both of those hold only while z = 0 lies inside the interval.  With z1 <= 0
    the whole domain sits in the bed: both branches of the `where` land below
    zero, the map folds back on itself at zeta = 0, and the result is [z0, 0]
    with [z1, 0] covered twice by reversed cells -- a domain that is not the one
    asked for, and that solves without complaining.  There is no interface to
    cluster towards in that case (nor for z0 >= 0, which folds the same way the
    other side up), so the mesh is left uniform instead.
    """
    domain = mesh.create_interval(MPI.COMM_SELF, n_cells, [0.0, 1.0])
    s = domain.geometry.x[:, 0]                            # uniform in [0, 1]

    if not z0 < 0.0 < z1:
        domain.geometry.x[:, 0] = z0 + (z1 - z0) * s
        return domain

    zeta = 2.0 * s - 1.0                                   # uniform in [-1, 1]
    f = (np.sign(zeta) * (np.exp(grading * np.abs(zeta)) - 1.0)
         / (np.exp(grading) - 1.0))
    domain.geometry.x[:, 0] = np.where(f >= 0.0, z1 * f, -z0 * f)
    return domain


# ===========================================================================
# Weak form and solve
# ===========================================================================
def solve_profile(domain, z0, z1, slope, d_p, eps_b, rho, nu, reg, g=G,
                  kappa=KAPPA, lam=LAM, A_star=A_STAR, A_E=A_E, B_E=B_E,
                  initial=None, options=None):
    """Solve the momentum balance on `domain`; return (zs, us), sorted by z.

    Multiplying the strong form by v and integrating tau' by parts gives

        F(u; v) = int_Omega [ source*v - tau*v' ] dz + [tau*v]_boundary = 0.

    The boundary term drops at both ends: at z1 because the free surface
    carries no stress, and at z0 because u is prescribed there (so v = 0).

    Newton starts from the uniform bed velocity u_SSL rather than from zero.
    At u = 0 both nonlinear stresses vanish -- van Driest leaves only
    ell = kappa*Z*sqrt(reg), and Forchheimer only sqrt(reg)*u -- so the first
    Jacobian is the *laminar* operator and its step aims at g*i*H^2/nu, which
    is 14 m/s for water against a true 0.4 m/s.  The line search cannot walk
    that back and reports DIVERGED_LINE_SEARCH on iteration 0.  Seeding with
    u_SSL switches the turbulence on in the very first Jacobian; the overshoot
    scales with 1/nu, which is why the zero guess survives nu = 3e-6 but not
    the 1e-6 of water.
    """
    V = fem.functionspace(domain, ("Lagrange", 1))
    z = ufl.SpatialCoordinate(domain)[0]
    u = fem.Function(V, name="u")
    v = ufl.TestFunction(V)
    du, dv = ufl.grad(u)[0], ufl.grad(v)[0]  # 1D: ufl.grad is a length-1 vector

    u_bed = u_subsurface(slope=slope, d_p=d_p, eps_b=eps_b, rho=rho, nu=nu,
                         g=g, A_E=A_E, B_E=B_E)
    if initial is None:
        u.x.array[:] = u_bed               # initial guess; see the docstring
    else:
        u.x.array[:] = initial(V.tabulate_dof_coordinates()[:, 0])

    src = source(z, u, du, slope=slope, d_p=d_p, eps_b=eps_b, rho=rho, nu=nu,
                 reg=reg, g=g, A_E=A_E, B_E=B_E, engine=ufl)
    tau = stresses(z, u, du, d_p=d_p, eps_b=eps_b, rho=rho, nu=nu, reg=reg,
                   kappa=kappa, lam=lam, A_star=A_star, engine=ufl)

    F = (src * v - sum(tau) * dv) * ufl.dx
    J = ufl.derivative(F, u)              # exact Newton Jacobian, from UFL

    bed = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], z0))
    bcs = [fem.dirichletbc(u_bed, bed, V)]

    options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",      # backtracking: the robustness that a
        "snes_linesearch_order": "2",      # fixed relaxation factor lacks
        # Tightening these below ~1e-9 does not buy accuracy: it asks the line
        # search to resolve a residual difference smaller than double precision
        # carries against an O(rho*g*i) balance, so it grinds inside roundoff
        # and reports DIVERGED_LINE_SEARCH on solves that have converged.  The
        # failures are knife-edge -- which of a sweep's points trip is not even
        # reproducible between runs -- so they read as physics and are not.
        "snes_rtol": 1.0e-9,
        "snes_stol": 1.0e-12,
        # Absolute floor at ~1e-10 of the driving force carried by the column.
        # It has to track rho*g*i: a fixed floor is either unreachable, in which
        # case the line search grinds around inside roundoff and reports
        # DIVERGED_LINE_SEARCH on an already converged solve, or loose enough to
        # accept garbage.
        "snes_atol": 1.0e-10 * rho * g * slope * (z1 - z0),
        "snes_max_it": 1000,
        "ksp_type": "preonly",
        "pc_type": "lu",
    } | (options or {})

    problem = NonlinearProblem(F, u, bcs=bcs, J=J,
                               petsc_options_prefix="vprof_",
                               petsc_options=options)
    problem.solve()

    reason = problem.solver.getConvergedReason()
    if reason <= 0:
        raise RuntimeError(f"SNES diverged (reason {reason}).")

    z_dofs = V.tabulate_dof_coordinates()[:, 0]
    order = np.argsort(z_dofs)
    return z_dofs[order], u.x.array[order]


# ===========================================================================
# Verification against the strong form, independently of the solver
# ===========================================================================
def momentum_error(zs, src, tau, z0, z1, slope, rho, g=G):
    """max|sum(tau) - int_z^z1 source| / (rho g i H).

    From 0 = source + tau' and tau(z1) = 0 at the free surface,

        tau(z) = int_z^z1 source dz' .

    The mismatch is dominated by the np.gradient estimate of u' on a P1 field,
    so it converges at first order under mesh refinement rather than reflecting
    the accuracy of the solve.
    """
    cumulative = np.concatenate(
        ([0.0], np.cumsum(0.5 * (src[1:] + src[:-1]) * np.diff(zs)))
    )
    tau_expected = cumulative[-1] - cumulative
    return np.abs(tau - tau_expected).max() / (rho * g * slope * (z1 - z0))


# ===========================================================================
# Plot
# ===========================================================================
def plot_profile(zs, us, dus, es, des, fs, tau_v, tau_t, tau_d, ze08,
                 z1, slope, rho, nu, g=G):
    """Four panels: profile, stresses, balance contributions, relative error."""

    fig, axes = plt.subplots(ncols=4, figsize=(11.3, 4), sharey=True)
    # fig.suptitle(
    #     "Numerical resolution of the steady-state, uniform flow velocity "
    #     "profile over a porous bed"
    #     "\n"
    #     r"$M_x' = \epsilon\varrho g i + (\tau_v+\tau_t+\tau_d)' + \epsilon f "
    #     r"- \varrho\nu\epsilon' u' = 0$"
    # )

    for i, ax in enumerate(axes):
        ax.axhline(0.0, ls="-", lw=0.5, c="k", alpha=0.5,
                   label="$z=0$" if i == 0 else None)
        ax.axhline(ze08, ls="--", lw=0.5, c="k", alpha=0.5,
                   label=r"$z(\epsilon=0.8)$" if i == 0 else None)
        ax.grid(alpha=0.3)

    # 1. the velocity profile, with the porosity on a twin axis behind it
    ax_eps = axes[0].twiny()
    ax_eps.set_zorder(axes[0].get_zorder() - 1)
    axes[0].patch.set_visible(False)
    l_eps, = ax_eps.plot(es, zs, "-.", c="C1", alpha=0.8,
                         label=r"Porosity $\epsilon$")
    ax_eps.set_xlabel(r"Porosity  $\epsilon$ [-]",
                      c=l_eps.get_color(), alpha=l_eps.get_alpha())
    ax_eps.dataLim.x0 = 0.0
    ax_eps.dataLim.x1 = 1.0
    [t.set_color(l_eps.get_color()) for t in ax_eps.xaxis.get_ticklabels()]

    l_u, = axes[0].plot(us, zs, label="Velocity $u$")
    m = zs > - ze08
    q = np.trapezoid(us[m]*es[m], zs[m])
    axes[0].set_xlabel("Velocity  $u$ [m/s]  "rf"($q = {q*1e3:.2g}~\mathrm{{l/s/m}}$)", c=l_u.get_color())
    axes[0].set_ylabel("Vertical coordinate  $z$ [m]")
    axes[0].dataLim.x0 = 0.0
    axes[0].dataLim.x1 = us.max()
    [t.set_color(l_u.get_color()) for t in axes[0].xaxis.get_ticklabels()]

    # one legend for both axes of the first panel
    lines, labels = axes[0].get_legend_handles_labels()
    lines_eps, labels_eps = ax_eps.get_legend_handles_labels()
    axes[0].legend(lines + lines_eps, labels + labels_eps)

    # 2. the stresses themselves
    axes[1].plot(tau_v + tau_t + tau_d, zs, c="k", label=r"$\sum \tau_{(\cdot)}$")
    l_v, = axes[1].plot(tau_v, zs, ":", label=r"$\tau_v$")
    l_t, = axes[1].plot(tau_t, zs, "-.", label=r"$\tau_t$")
    l_d, = axes[1].plot(tau_d, zs, "--", label=r"$\tau_d$")
    axes[1].set_xlabel(r"Shear stress  $\tau$ [Pa]")
    axes[1].legend()

    def style(line):
        return {"ls": line.get_linestyle(), "c": line.get_color(),
                "alpha": line.get_alpha()}

    # 3. the objective, term by term: 0 = source + tau'
    axes[2].plot(np.gradient(tau_v, zs), zs, label=r"$\tau_v'$", **style(l_v))
    axes[2].plot(np.gradient(tau_t, zs), zs, label=r"$\tau_t'$", **style(l_t))
    axes[2].plot(np.gradient(tau_d, zs), zs, label=r"$\tau_d'$", **style(l_d))
    axes[2].plot(es * fs, zs, label=r"$\epsilon f$", lw=1)
    axes[2].plot(es * rho * g * slope, zs, c="k", label=r"$\epsilon \varrho g i$")
    axes[2].set_xlabel("Momentum balance contributions  [Pa/m]")
    axes[2].legend()

    # 4. the same residual, scaled by the driving force so that 1 = 100% off
    tau = tau_v + tau_t + tau_d
    z_lim = z1 - 1e-3
    m = zs < z_lim
    rerr = 100 * ((np.gradient(tau, zs) + es * fs - rho * nu * des * dus)
                  / (es * rho * g * slope) + 1)
    axes[3].plot(rerr[m], zs[m])
    axes[3].set_xlim(axes[3].get_xlim())
    axes[3].plot(rerr[~m], zs[~m], ls="--")
    axes[3].set_xlabel(r"Relative error  $\frac{\boldsymbol{p}_x'}{\epsilon\varrho g i}$ [%]")

    fig.tight_layout()
    return fig, axes


# ===========================================================================
# The same thing, with the parameters carried instead of passed
# ===========================================================================
@dataclass
class VelocityProfileSolver:
    """One run: every function above, bound to a single set of parameters.

    The functions stay the reference -- each method here does nothing but call
    one of them with this case's numbers.  The defaults are the UPPERCASE
    block, so ``Case()`` is the reference case and ``Case(slope=3e-2)`` is that
    case with one number changed::

        zs, us = Case(slope=3.0e-2).run()

    The two quantities that follow from the parameters alone -- the interface
    sharpness a and the deep-bed velocity u_bed -- are computed once in
    __post_init__ and are not constructor arguments.
    """

    # Solver settings
    n_cells: int = N_CELLS
    grading: float = GRADING
    reg: float = REG

    # Physical constants
    g: float = G
    kappa: float = KAPPA
    lam: float = LAM
    A_star: float = A_STAR
    A_E: float = A_E
    B_E: float = B_E

    def update(self, other_dict):
        self.__dict__.update(other_dict)

    @property
    def a(self):
        return interface_sharpness(self.d_p)

    @property
    def u_bed(self):
        return u_subsurface(
            slope=self.slope, d_p=self.d_p, eps_b=self.eps_b, rho=self.rho,
            nu=self.nu, g=self.g, A_E=self.A_E, B_E=self.B_E
        )

    def z_from_e(self, e):
        return 2 / self.a * np.arctanh(2 * (e - self.eps_b) / (1 - self.eps_b) - 1)

    @property
    def ze08(self):
        return self.z_from_e(0.8)

    @property
    def depth(self):
        return self.z1 - self.ze08

    # -- profiles ----------------------------------------------------------
    def porosity(self, z, engine=np):
        return porosity(z, d_p=self.d_p, eps_b=self.eps_b, engine=engine)

    def porosity_rate(self, z, engine=np):
        return porosity_rate(z, d_p=self.d_p, eps_b=self.eps_b, engine=engine)

    def integral_depth(self, z, engine=np):
        return integral_depth(z, d_p=self.d_p, engine=engine)

    def mixing_length(self, z, u, engine=np):
        return mixing_length(z, u, d_p=self.d_p, nu=self.nu, reg=self.reg,
                             kappa=self.kappa, A_star=self.A_star,
                             engine=engine)

    # -- stresses and body forces ------------------------------------------
    def stresses(self, z, u, du, engine=np):
        return stresses(z, u, du, d_p=self.d_p, eps_b=self.eps_b, rho=self.rho,
                        nu=self.nu, reg=self.reg, kappa=self.kappa,
                        lam=self.lam, A_star=self.A_star, engine=engine)

    def drag(self, e, u, engine=np):
        return drag(e, u, d_p=self.d_p, rho=self.rho, nu=self.nu, reg=self.reg,
                    A_E=self.A_E, B_E=self.B_E, engine=engine)

    def source(self, z, u, du, engine=np):
        return source(z, u, du, slope=self.slope, d_p=self.d_p,
                      eps_b=self.eps_b, rho=self.rho, nu=self.nu, reg=self.reg,
                      g=self.g, A_E=self.A_E, B_E=self.B_E, engine=engine)

    # -- mesh and solve ----------------------------------------------------
    def mesh(self):
        """A fresh graded interval; `solve` builds one if given none."""
        return graded_mesh(z0=self.z0, z1=self.z1, n_cells=self.n_cells,
                           grading=self.grading)

    def solve_naive(self, domain=None):
        """(zs, us), sorted by z.  Pass `domain` to reuse an existing mesh."""
        if domain is None:
            domain = self.mesh()
        return solve_profile(domain, z0=self.z0, z1=self.z1, slope=self.slope,
                             d_p=self.d_p, eps_b=self.eps_b, rho=self.rho,
                             nu=self.nu, reg=self.reg, g=self.g,
                             kappa=self.kappa, lam=self.lam,
                             A_star=self.A_star, A_E=self.A_E, B_E=self.B_E)

    def solve_initial(self, domain=None):
        """(zs, us), sorted by z.  Pass `domain` to reuse an existing mesh."""
        return solve_profile(domain, z0=self.z0, z1=self.z1, slope=self.slope,
                             d_p=self.d_p, eps_b=self.eps_b, rho=self.rho,
                             nu=self.nu, reg=self.reg, g=self.g,
                             kappa=0., lam=0.,
                             A_star=self.A_star, A_E=self.A_E, B_E=self.B_E)

    def solve(self, domain=None):
        """(zs, us), sorted by z.  Pass `domain` to reuse an existing mesh."""
        if domain is None:
            domain = self.mesh()
        zs_initial, us_initial = self.solve_initial(domain)
        return solve_profile(domain, initial=lambda z: np.interp(z, zs_initial, us_initial),
                             z0=self.z0, z1=self.z1, slope=self.slope,
                             d_p=self.d_p, eps_b=self.eps_b, rho=self.rho,
                             nu=self.nu, reg=self.reg, g=self.g,
                             kappa=self.kappa, lam=self.lam,
                             A_star=self.A_star, A_E=self.A_E, B_E=self.B_E)

    # -- post-processing ---------------------------------------------------
    def fields(self, zs, us):
        """Every derived field the checks and the plot need, as a dict."""
        dus = np.gradient(us, zs)
        es = self.porosity(zs)
        des = self.porosity_rate(zs)
        fs = self.drag(es, us)
        tau_v, tau_t, tau_d = self.stresses(zs, us, dus)
        return {
            "dus": dus, "es": es, "des": des, "fs": fs,
            "tau_v": tau_v, "tau_t": tau_t, "tau_d": tau_d,
            "src": self.source(zs, us, dus),
            "ell": self.mixing_length(zs, us),
            "Z": self.integral_depth(zs),
            "ze08": self.ze08,
        }

    def momentum_error(self, zs, src, tau):
        return momentum_error(zs, src, tau, z0=self.z0, z1=self.z1,
                              slope=self.slope, rho=self.rho, g=self.g)

    def plot(self, zs, us, f=None):
        """Four panels; `f` reuses a `fields` dict instead of recomputing."""
        f = self.fields(zs, us) if f is None else f
        return plot_profile(zs, us, f["dus"], f["es"], f["des"], f["fs"],
                            f["tau_v"], f["tau_t"], f["tau_d"], f["ze08"],
                            z1=self.z1, slope=self.slope, rho=self.rho,
                            nu=self.nu, g=self.g)

    # -- driver ------------------------------------------------------------
    def describe(self, domain):
        """What is about to be solved -- printed before the solve, not after."""
        dz_cells = np.diff(np.sort(domain.geometry.x[:, 0]))
        print(f"u_SSL(i={self.slope:g}) = {self.u_bed * 1e3:.3f} mm/s")
        print(f"mesh: {self.n_cells} cells, dz in "
              f"[{dz_cells.min() * 1e3:.4f}, {dz_cells.max() * 1e3:.4f}] "
              f"mm, interface width 1/a = {1e3 / self.a:.3f} mm")

    def check(self, zs, us, f):
        """The three diagnostics of a finished solve; returns the error."""
        tau = f["tau_v"] + f["tau_t"] + f["tau_d"]
        error = self.momentum_error(zs, f["src"], tau)
        h = self.z1 - f["ze08"]
        um = us[zs >= f["ze08"]].mean()
        print("momentum check: max|sum(tau) - expected| / (rho g i H) = "
              f"{error:.3e}")
        print(f"u in [{us.min():.4e}, {us.max():.4e}] m/s, "
              f"q = {np.trapezoid(us, zs) * 1e3:.4f} l/s/m")
        print(f"um / sqrt(g h i) = {um / np.sqrt(self.g * h * self.slope):.2f}"
              f" | h / d_p = {h / self.d_p:.2f}")
        return error

    def run(self, verbose=False):
        """Solve, check, optionally plot; return (zs, us)."""
        domain = self.mesh()
        if verbose:
            self.describe(domain)

        zs, us = self.solve(domain)
        f = self.fields(zs, us)
        if verbose:
            self.check(zs, us, f)

        return zs, us


# ===========================================================================
# Driver
# ===========================================================================
def main():
    """The reference case: the UPPERCASE block is `Case`'s set of defaults."""
    solver = VelocityProfileSolver()
    solver.update({
        # Variables of the case
        "rho": RHO,
        "nu": NU,
        "slope": SLOPE,
        "d_p": D_P,
        "eps_b": EPS_B,
        "z0": Z0,
        "z1": Z1,
    })
    zs, us = solver.run()
    solver.plot(zs, us, solver.fields(zs, us))
    plt.show()
    plt.plot(us, zs)
    plt.box(False)
    plt.savefig("standalone.svg")
    plt.show()


if __name__ == "__main__":
    main()
