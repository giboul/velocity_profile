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

Every number of a run -- physical constants, variables of the case, solver
settings -- lives in one `Parameters` object, `p`, which each function takes
instead of a list of individual arguments.  Its defaults are the UPPERCASE
block below, so `Parameters()` is the reference case and
`Parameters(slope=3e-2)` is that case with one number changed.
`p.solve()` returns a `Result`: the same parameters plus the solution and
every field derived from it.

The profiles also take `engine`, the module their maths goes through.  It
defaults to `numpy`, so calling them plainly evaluates them on arrays; the weak
form is the one place that passes `ufl` instead.

Run:  python momentum_ODE_solver.py
"""

from dataclasses import dataclass, replace
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ufl
from dolfinx import fem
from dolfinx import mesh as dx_mesh
from dolfinx.fem.petsc import NonlinearProblem
from mpi4py import MPI

style = Path(__file__).parent / "communications/hyporheic_flow/lhe.mplstyle"
plt.style.use(style)
# Consistent naming between the two engines, so `engine.<f>` works for both.
ufl.maximum = ufl.max_value
ufl.log = ufl.ln


# ===========================================================================
# Parameters -- the reference case, used as the defaults of `Parameters`
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
Z1 = 0.03               # free surface                         [m]
Z0 = -12e-3            # bottom of the bed                    [m]

# Solver settings
N_CELLS = 10_000       # 1D, so we can afford to be generous
REFINEMENT = 10.0      # cells at the interface are this much smaller than far
# from it (1 = uniform); see `refined_mesh`
REG = 1.0e-12          # smooths |.| and sqrt(.) at the origin so that their
# derivatives stay finite for Newton


@dataclass
class Parameters:
    """Every number of one run; the `p` that all the functions below take.

    The fields are plain attributes, so a sweep can change one in place
    (``p.z1 = ...``), and ``dataclasses.replace(p, kappa=0.0)`` gives a
    modified copy.  The quantities that follow from the parameters alone are
    properties, so they always match the current values.
    """

    # Variables of the case
    rho: float = RHO
    nu: float = NU
    slope: float = SLOPE
    d_p: float = D_P
    eps_b: float = EPS_B
    z0: float = Z0
    z1: float = Z1

    # Physical constants
    g: float = G
    kappa: float = KAPPA
    lam: float = LAM
    A_star: float = A_STAR
    A_E: float = A_E
    B_E: float = B_E

    # Solver settings
    n_cells: int = N_CELLS
    refinement: float = REFINEMENT
    reg: float = REG

    # -- derived quantities --------------------------------------------------
    @property
    def a(self):
        """a [1/m]: eps goes eps_b -> 99.9% of 1 over roughly one grain diameter."""
        return -2.0 / self.d_p * np.log(1.0 / 0.994 - 1.0)

    @property
    def u_bed(self):
        """Deep-bed equilibrium: the positive root of eps*rho*g*i + eps*fbar = 0."""
        e = self.eps_b
        qa = - self.rho * self.B_E * e * (1.0 - e) / self.d_p
        qb = - self.rho * self.nu * self.A_E * (1.0 - e)**2 / (e * self.d_p**2)
        qc = + e * self.rho * self.g * self.slope
        return (-qb - np.sqrt(qb**2 - 4.0 * qa * qc)) / (2.0 * qa)

    def z_from_e(self, e):
        """The height at which the porosity equals `e`."""
        return 2 / self.a * np.arctanh(2 * (e - self.eps_b) / (1 - self.eps_b) - 1)

    @property
    def ze08(self):
        """The height at which eps = 0.8: the conventional top of the bed."""
        return self.z_from_e(0.8)

    @property
    def depth(self):
        """Water depth h, from z(eps = 0.8) to the free surface."""
        return self.z1 - self.ze08

    # -- solving -------------------------------------------------------------
    def solve(self):
        """Solve this case (see `solve`) and return it as a `Result`."""
        return Result(**self.__dict__, **solve(self, mesh(self)))


# ===========================================================================
# Vertical profiles
# ===========================================================================


def porosity(z, p, engine=np):
    """eps(z), a smooth step from eps_b in the bed to 1 in the free water."""
    return p.eps_b + (1.0 - p.eps_b) * 0.5 * (1.0 + engine.tanh(0.5 * p.a * z))


def porosity_rate(z, p, engine=np):
    """eps'(z), analytically: (1-eps_b)*a/4 * sech^2(a*z/2)."""
    return (1.0 - p.eps_b) * p.a / 4.0 * (1.0 - engine.tanh(0.5 * p.a * z)**2)


def integral_depth(z, p, engine=np):
    """Z(z) = softplus(a*z)/a: height above the interface, smoothed over 1/a.

    Written as max(x, 0) + log(1 + exp(-|x|)) rather than log(1 + exp(x)):
    the two are the same function, but the second overflows float64 at
    a*z = 709, i.e. z = 51.3*d_p.  The inf then reaches tau_t as inf - inf
    and SNES aborts on iteration 0 with a NaN residual (reason -4).
    """
    return engine.maximum(z, 0.0) + engine.log(1.0 + engine.exp(-abs(p.a * z))) / p.a


def mixing_length(z, u, p, engine=np):
    """Prandtl length kappa*Z with van Driest near-wall damping."""
    Z = integral_depth(z, p, engine=engine)
    u_pos = engine.maximum(u, 0.0)
    damping = 1.0 - \
        engine.exp(-engine.sqrt(Z * u_pos / (p.nu * p.A_star**2) + p.reg))
    return p.kappa * Z * damping


# ===========================================================================
# Stresses and body forces
# ===========================================================================


def stresses(z, u, du, p, engine=np):
    """(tau_v, tau_t, tau_d) -- viscous, turbulent and dispersive."""
    eps = porosity(z, p, engine=engine)
    deps = porosity_rate(z, p, engine=engine)
    ell = mixing_length(z, u, p, engine=engine)
    tau_v = p.rho * p.nu * (deps * u + eps * du)
    tau_t = p.rho * eps * ell**2 * engine.sqrt(du * du + p.reg) * du
    tau_d = p.rho * p.lam * p.d_p / \
        (1.0 - p.eps_b) * (1.0 - eps) * eps * u * du
    return tau_v, tau_t, tau_d


def drag(e, u, p, engine=np):
    """Ergun resistance eps*fbar: Kozeny-Carman plus Forchheimer."""
    return (
        - p.A_E * (1.0 - e)**2 / (e**2 * p.d_p**2) * p.rho * p.nu * u
        - p.B_E * (1.0 - e) / p.d_p * p.rho * engine.sqrt(u * u + p.reg) * u
    )


def source(z, u, du, p, engine=np):
    """Every term of the strong form except tau', i.e. 0 = source(u) + tau'."""
    e = porosity(z, p, engine=engine)
    deps = porosity_rate(z, p, engine=engine)
    f = drag(e, u, p, engine=engine)
    return e * p.rho * p.g * p.slope + e * f - p.rho * p.nu * deps * du


# ===========================================================================
# Mesh
# ===========================================================================


def uniform_mesh(p):
    return dx_mesh.create_interval(MPI.COMM_SELF, p.n_cells, [p.z0, p.z1])


def refined_mesh(p):
    """A mesh that follows the porosity rate: the cell size is ~ 1/w(z) with

        w(z) = 1 + (refinement - 1) * eps'(z) / eps'(0),

    so the cells at the interface are `p.refinement` times smaller than far
    from it.  eps' is taken over twice the grain diameter so that the
    refinement covers the whole transition.  The nodes sit at equal
    increments of W = int w, found by interpolation on a fine sample.
    """
    wide = replace(p, d_p=2 * p.d_p)
    z = np.linspace(p.z0, p.z1, 50 * p.n_cells + 1)
    w = 1 + (p.refinement - 1) * porosity_rate(z, wide) / \
        porosity_rate(0, wide)
    W = np.cumsum(w) * (z[1] - z[0])
    nodes = np.interp(np.linspace(W[0], W[-1], p.n_cells + 1), W, z)

    # map the uniform nodes onto the graded ones, whatever order dolfinx
    # stores the geometry in
    domain = uniform_mesh(p)
    x = domain.geometry.x
    x[:, 0] = np.interp(x[:, 0], np.linspace(p.z0, p.z1, p.n_cells + 1), nodes)
    return domain


mesh = refined_mesh

# ===========================================================================
# Weak form and solve
# ===========================================================================


def solve_profile(domain, p, initial=None, options=None):
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
    du, dv = ufl.grad(u)[0], ufl.grad(v)[0]    # 1D: ufl.grad is a 1-vector

    u_bed = p.u_bed
    if initial is None:
        u.x.array[:] = u_bed                   # initial guess; see docstring
    else:
        u.x.array[:] = initial(V.tabulate_dof_coordinates()[:, 0])

    src = source(z, u, du, p, engine=ufl)
    tau = stresses(z, u, du, p, engine=ufl)

    F = (src * v - sum(tau) * dv) * ufl.dx
    # exact Newton Jacobian, from UFL
    J = ufl.derivative(F, u)

    bed = fem.locate_dofs_geometrical(V, lambda x: np.isclose(x[0], p.z0))
    bcs = [fem.dirichletbc(u_bed, bed, V)]

    options = {
        "snes_type": "newtonls",
        "snes_linesearch_type": "bt",
        "snes_linesearch_order": "2",
        "snes_rtol": 1.0e-9,
        "snes_stol": 1.0e-12,
        "snes_atol": 1.0e-10 * p.rho * p.g * p.slope * (p.z1 - p.z0),
        "snes_max_it": 1000,
        "ksp_type": "preonly",
        "pc_type": "lu",
    } | (options or {})

    problem = NonlinearProblem(
        F, u, bcs=bcs, J=J, petsc_options=options, petsc_options_prefix="vprof_")
    problem.solve()

    reason = problem.solver.getConvergedReason()
    if reason <= 0:
        raise RuntimeError(f"SNES diverged (reason {reason}).")

    z_dofs = V.tabulate_dof_coordinates()[:, 0]
    order = np.argsort(z_dofs)
    return z_dofs[order], u.x.array[order]


def solve_naive(p, domain=None):
    """(zs, us), sorted by z, in one Newton run from u_SSL.

    Pass `domain` to reuse an existing mesh.
    """
    if domain is None:
        domain = mesh(p)
    return solve_profile(domain, p)


def solve_initial(p, domain):
    """(zs, us) without turbulence or dispersion: the initial guess of `solve`."""
    return solve_profile(domain, replace(p, kappa=0.0, lam=0.0))


def solve(p, domain=None):
    """(z, u), sorted by z, continued from the laminar solution.

    Pass `domain` to reuse an existing mesh.
    """
    if domain is None:
        domain = mesh(p)
    z_initial, u_initial = solve_initial(p, domain)
    z, u = solve_profile(
        domain, p, initial=lambda z: np.interp(z, z_initial, u_initial))
    return {"z": z, "u": u}


# ===========================================================================
# Post-processing
# ===========================================================================


@dataclass
class Result(Parameters):
    """A solved profile: the parameters it was solved with, plus every field
    derived from the solution (z, u).

    Being a `Parameters`, it can be passed wherever a `p` is expected, so the
    profile functions above work on it directly.  Build one with
    ``p.solve()``.  The derived arrays are computed once, in
    ``__post_init__``, and are plain attributes afterwards; z and u are
    stored sorted by z.
    """

    z: np.ndarray = None                       # nodes, sorted           [m]
    u: np.ndarray = None                       # velocity at the nodes   [m/s]

    def __post_init__(self):
        order = np.argsort(self.z)
        self.z = np.asarray(self.z, dtype=float)[order]
        self.u = np.asarray(self.u, dtype=float)[order]

        self.du = np.gradient(self.u, self.z)              # u'           [1/s]
        self.e = porosity(self.z, self)                    # eps          [-]
        self.de = porosity_rate(self.z, self)              # eps'         [1/m]
        # fbar         [Pa/m]
        self.f = drag(self.e, self.u, self)
        self.tau_v, self.tau_t, self.tau_d = \
            stresses(self.z, self.u, self.du, self)        # stresses     [Pa]
        self.tau = self.tau_v + self.tau_t + self.tau_d
        # tau'         [Pa/m]
        self.dtau = np.gradient(self.tau, self.z)
        self.src = source(self.z, self.u, self.du, self)   # 0 = src + tau'
        self.ell = mixing_length(self.z, self.u, self)     # mixing length [m]
        self.Z = integral_depth(self.z, self)              # smoothed height

        water = self.z > -self.ze08                        # above the bed
        self.u_mean = self.u[water].mean()                 # [m/s]
        self.q = np.trapezoid(self.u[water] * self.e[water], self.z[water])
        self.froude = self.u_mean / np.sqrt(self.g * self.depth * self.slope)

    def momentum_error(self, relative=True):
        return momentum_error(self, relative=relative)

    def total_momentum_error(self, relative=True):
        return np.sum(self.momentum_error(relative=relative) * np.gradient(self.z))

    def plot(self):
        return plot_profile(self)


def momentum_error(r, relative=True):
    """Pointwise residual of the strong form, relative to eps*rho*g*i.

    Built from nodal finite differences of a P1 field, so it is only an
    indication: it is first order at the end nodes, meaningless at the
    Dirichlet node, and O(1) over the last few nodes below the free
    surface, where u' ~ sqrt(z1 - z) makes u'' unbounded.
    """
    driving = r.e * r.rho * r.g * r.slope
    heterogeneity = r.rho * r.nu * r.de * r.du
    err = driving + r.dtau + r.e * r.f - heterogeneity
    return err / driving if relative else err


# ===========================================================================
# Plot
# ===========================================================================


def plot_profile(r):
    """Four panels: profile, stresses, balance contributions, relative error.

    `r` is a `Result`.
    """
    zs, us, dus, es, des, fs = r.z, r.u, r.du, r.e, r.de, r.f
    tau_v, tau_t, tau_d, ze08 = r.tau_v, r.tau_t, r.tau_d, r.ze08
    slope, rho, nu, g = r.slope, r.rho, r.nu, r.g

    fig, axes = plt.subplots(ncols=4, figsize=(11.3, 4), sharey=True,
                             gridspec_kw={"wspace": 0})

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
    for t in ax_eps.xaxis.get_ticklabels():
        t.set_color(l_eps.get_color())

    l_u, = axes[0].plot(us, zs, label="Velocity $u$")
    axes[0].set_xlabel(f"Velocity  $u$ [m/s], $q = {r.q * 1e3:.2g}$ [l/s/m]",
                       c=l_u.get_color())
    axes[0].set_ylabel("Vertical coordinate  $z$ [m]")
    axes[0].dataLim.x0 = 0.0
    axes[0].dataLim.x1 = us.max()
    for t in axes[0].xaxis.get_ticklabels():
        t.set_color(l_u.get_color())

    # one legend for both axes of the first panel
    lines, labels = axes[0].get_legend_handles_labels()
    lines_eps, labels_eps = ax_eps.get_legend_handles_labels()
    axes[0].legend(lines + lines_eps, labels + labels_eps)

    # 2. the stresses themselves
    axes[1].plot(r.tau, zs, c="k", label=r"$\sum \tau_{(\cdot)}$")
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
    axes[2].plot(rho * nu * (des * us + es * dus), zs,
                 label=r"$(\epsilon u)'\mu$", lw=1)
    axes[2].plot(es * fs, zs, label=r"$\epsilon f$", lw=1)
    axes[2].plot(es * rho * g * slope, zs, c="k",
                 label=r"$\epsilon \varrho g i$")
    axes[2].set_xlabel("Momentum balance contributions  [Pa/m]")
    axes[2].legend()

    # 4. the pointwise residual of the strong form, scaled by the driving
    #    force so that 100 = 100% off (see `Result.momentum_error`)

    ax_err = axes[3].twiny()
    ax_err.xaxis.tick_bottom()
    ax_err.xaxis.set_label_position("bottom")
    axes[3].xaxis.tick_top()
    axes[3].xaxis.set_label_position("top")

    fill_color = axes[3].fill_betweenx(
        zs, np.gradient(zs), alpha=0.3).get_fc()
    fill_color[:, -1] = 1.
    axes[3].set_xlabel("Cell size [m]", c=fill_color)
    axes[3].set_xlim(0, None)
    for t in axes[3].xaxis.get_ticklabels():
        t.set_color(fill_color)

    m = 4
    ax_err.plot(100 * r.momentum_error()[m:-m], zs[m:-m], "-k", mfc="none")
    ax_err.set_xlabel(
        r"Relative residual  $\frac{\boldsymbol{p}_x'}{\epsilon\varrho g i}$ [%]")

    axes[0].set_ylim(r.z0, r.z1)

    fig.tight_layout()
    return fig, axes


# ===========================================================================
# Driver
# ===========================================================================


def main():
    r = Parameters().solve()
    print(f" | e = {r.total_momentum_error():.3g} %"
          f" | q = {r.q * 1e3:.4f} l/s/m"
          f" | u_mean / sqrt(g h i) = {r.froude:.2f}"
          f" | h / d_p = {r.depth / r.d_p:.2f}"
          f" | max |residual| = {np.abs(r.momentum_error()).max():.3e}")
    r.plot()
    print()
    plt.show()


if __name__ == "__main__":
    main()
