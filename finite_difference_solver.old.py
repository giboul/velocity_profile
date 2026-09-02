r"""
Solving the following equation for a particular porosity function

$$0 = \epsilon\varrho g i + \epsilon\bar f + \left(\tau_v + \tau_t + \tau_d\right)' - \varrho\nu\epsilon'u'$$
$$\tau_v = \varrho\nu\left(\epsilon u\right)'$$
$$\tau_t = -\varrho\epsilon\kappa^2 Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2$$
$$\tau_d = - \frac{\varrho\lambda d_p}{1-\epsilon_b}[1-\epsilon]\epsilon uu'$$
"""
import numpy as np
from matplotlib import pyplot as plt
from numpy.polynomial import Polynomial
from scipy.linalg import solve_banded

# np.seterr(all="raise")


g = 9.81  # kg m / s²
kappa = 0.41  # [-]
lam = 0.015  # [-]

rho = 950.0  # kg / m³
nu = 3.1e-6  # m²/s

i = 1.0 / 100
dp = 8. / 1000
eps_bulk = 0.5

a = -2/dp * np.log(1/(99/100) - 1)
print(f"{a = }")

A_E = 180.0
B_E = 1.75
A = 26.0

start = - 0.5
end = + 0.3
n = 100_000

# Porosity functions

def porosity(z):
    return eps_bulk + (1.0 - eps_bulk) / (1.0 + np.exp(-a * z))

def no_overflow_reduc(z):

    if isinstance(z, float):
        cast_float = True
        z = np.array([z])
    else:
        cast_float = False

    exp = np.where(0 < z, np.exp(-a * z), np.exp(a * z))
    r = exp / (1.0 + exp)**2

    if cast_float is True:
        return r[0]
    return r

def porosity_rate(z):
    r = no_overflow_reduc(z)
    out = (1.0 - eps_bulk) * a * r
    return out

def porosity_accel(z):
    exp = np.exp(-a * z)
    r = no_overflow_reduc(z)
    return a**2 * (1 - eps_bulk) * r * (2 * exp / (1.0 + exp) - 1)

def integral_depth(z):
    return np.log1p(np.exp(a * z)) / a

def integral_depth_rate(z):
    return np.exp(a * z) / (1.0 + np.exp(a * z))

# Subsurface flow

def kozeny_carman(u, e):
    return (
        - rho * nu * A_E * (1 - e)**2 / (e**2 * dp**2) * u
        - rho * B_E * (1 - e) / dp * u**2
    )

def steady_state_subsurface_flow_speed(e):
    a = - rho * B_E * e * (1 - e) / dp
    b = - rho * nu * A_E * (1 - e)**2 / (e * dp**2)
    c = + e * rho * g * i
    # print(f"{Polynomial([c, b, a]).roots().max() = :.2e}")
    return (- b - np.sqrt(b**2 - 4 * a * c)) / (2 * a)

u_SSL = steady_state_subsurface_flow_speed(eps_bulk)
K_KC_Darcy = dp**2 * eps_bulk**3 / (A_E * (1-eps_bulk)**2)
print(f"u_SSL = {u_SSL*1e3:.2f} mm/s ~ {rho*g*i/(rho*nu)*K_KC_Darcy*1e3:.2f} mm/s")
print(f"K_KC  = {rho*nu*eps_bulk/(rho*g*i)*u_SSL*1e6:.3f} mm² ~ {K_KC_Darcy*1e6:.3f} mm²")

# ODE

## Factors in the ODE

def d0u0_factor(e):
    return e*rho*g*i

def d0u1_factor(e, d2e):
    return rho*nu*d2e - rho*nu*A_E*(1-e)**2/(e*dp**2)

def d0u2_factor(e):
    return - rho * B_E * (1 - e) * e / dp

def d0u_factor(u, e, d2e):
    return d0u0_factor(e) + d0u1_factor(e, d2e)*u + d0u2_factor(e)*u**2

def d1u_factor(u, e, de):
    return rho*nu*de + rho*dp*lam/(1-eps_bulk)*(1-2*e)*de*u

def d1u2_factor(u, zls, dzls, e, de):
    phi = np.sqrt(zls * u / (nu * A**2))
    exp = np.exp(-phi)
    return + rho * kappa**2 * (
        + de * zls**2 * (1 - exp)**2
        + 2 * e * zls * dzls * (1 - exp)**2
        + e * zls**2 * (1 - exp) * exp * dzls * u / np.sqrt(nu * A**2 * zls * u)
    ) + rho * dp * lam / (1 - eps_bulk) * (1 - e) * e

def d1u3_factor(u, zls, e):
    phi = np.sqrt(zls * u / (nu * A**2))
    exp = np.exp(-phi)
    return + rho * kappa**2 * e * zls**2 * (1 - exp) * exp * zls / np.sqrt(nu * A**2 * zls * u)

def d2u_factor(u, e):
    return rho * nu * e - rho * dp * lam / (1 - eps_bulk) * (1 - e) * e * u

def d2u_d1u_factor(u, zls, e):
    exp = np.exp(-np.sqrt(zls * u / (nu * A**2)))
    return + 2 * rho * kappa**2 * e * zls**2 * (1 - exp)**2

## Derivatives computers

def du_computer(z, u, d2u, debug=False):

    e = porosity(z)
    de = porosity_rate(z)
    d2e = porosity_accel(z)
    zls = integral_depth(z)
    dzls = integral_depth_rate(z)

    coefs = np.array([
        d0u_factor(u, e, d2e) + d2u_factor(u, e) * d2u,
        d1u_factor(u, e, de) + d2u_d1u_factor(u, zls, e) * d2u,
        d1u2_factor(u, zls, dzls, e, de),
        d1u3_factor(u, zls, e),
    ])

    roots = Polynomial(coefs).roots()
    root = roots.real.max()

    if debug is True:
        print("="*20)
        print("0 = " + " ".join([f"{c:+.2e}{xn}" for xn, c in zip(('', '*x', '*x²', '*x³'), coefs)]))
        print(f"{roots = }")
        zero_check = (coefs[:, None] * roots[None, :] ** np.arange(coefs.size)[:, None]).sum(axis=0)
        print("f(roots) = " + " | ".join([f"{z:.2e}" for z in zero_check]))
        print(f"RR : {root = :.2e}")

    return roots.real.max()

def d2u_computer(z, u, du, debug=False):

    e = porosity(z)
    de = porosity_rate(z)
    d2e = porosity_accel(z)
    zls = integral_depth(z)
    dzls = integral_depth_rate(z)

    p = d2u_factor(u, e) + d2u_d1u_factor(u, zls, e) * du

    q = (
        + d0u_factor(u, e, d2e)
        + d1u_factor(u, e, de) * du
        + d1u2_factor(u, zls, dzls, e, de) * du**2
        + d1u3_factor(u, zls, e) * du**3
    )

    if debug is True:
        print("#"*20)
        print(f"d2u       <= {d2u_factor(u, e):+.2e}")
        print(f"d2u * d1u <= {d2u_d1u_factor(u, zls, e):+.2e}")
        print(f"d1u^3     <= {d1u3_factor(u, zls, e):+.2e}")
        print(f"d1u^2     <= {d1u2_factor(u, zls, dzls, e, de):+.2e}")
        print(f"d1u       <= {d1u_factor(u, e, de):+.2e}")
        print(f"d0u       <= {d0u_factor(u, e, d2e):+.2e}")
        print(f"=> {- q / p = :.2e}")

    return - q / p

def friedli(z, u0=u_SSL, uf=None, bound_u0=False, debug=False, live_fig=True):

    u = np.full(z.shape, u0, dtype=np.float64)
    du = 0
    e = porosity(z)
    de = porosity_rate(z)
    d2e = porosity_accel(z)
    zls = integral_depth(z)
    dzls = integral_depth_rate(z)

    main_band = np.zeros(z.size)
    upper_band = np.zeros(z.size)
    lower_band = np.zeros(z.size)
    grad_band = np.zeros(z.size)

    # upper_band[0] = np.nan
    # lower_band[-1] = np.nan
    # grad_band[-2:] = np.nan

    nb = 2
    nt = 2

    if live_fig:
        fig, (ax, ax_e) = plt.subplots(ncols=2, sharey=True)
        line, = ax.plot(u, z, '-_')
        ax.axline((u_SSL, start), slope=np.inf, ls="-.", alpha=0.5, c="k")
        ax.set_xlabel("$u$", c="C0")
        ax.set_ylabel("$z$")
        ax.dataLim.x0 = 0.
        ax.dataLim.x1 = 5*u0

        # ax_e = ax.twiny()
        ax_e.plot(porosity(z), z, '-.', c="C1", alpha=0.5)
        ax_e.axline((eps_bulk, start), slope=np.inf, ls=":", c="k", alpha=0.5)
        ax_e.axline((1, start), slope=np.inf, ls=":", c="k", alpha=0.5)
        ax_e.axline((0, 0), slope=0, ls=":", c="k", alpha=0.5)
        ax_e.set_xlabel(r"$\epsilon$", c="C1")

        ax_bis = ax_e.twiny()
        line_bis, = ax_bis.plot(u, z, '-_')

        fig.show(warn=True)

    for iteration in range(10):
        C = - d0u0_factor(e)
        L = d0u1_factor(e, d2e) + d0u2_factor(e)*u
        S = d1u_factor(u, e, de) + d1u2_factor(u, zls, dzls, e, de)*du + d1u3_factor(u, zls, e)*du**2
        A = d2u_factor(u, e) + d2u_d1u_factor(u, zls, e) * du

        print()
        print(f"{C[1] = :10.2e} | {C[-5] = :10.2e}")
        print(f"{L[1] = :10.2e} | {L[-5] = :10.2e}")
        print(f"{S[1] / (2 * dz) = :10.2e} | {S[-5] / (2 * dz) = :10.2e}")
        print(f"{A[1]/(4* dz**2) = :10.2e} | {A[-5]/(4* dz**2) = :10.2e}")
        print()

        # Assembling matrix

        ## Governing equation

        ## L * u_{0}
        main_band[:] = L

        ### S * (u_{1} - u_{-1}) / (2 * Δz)
        upper_band[+1:] = - S[:-1] / (2 * dz)
        lower_band[:-1] = + S[+1:] / (2 * dz)

        upper_band[+1:] += + A[:-1] / dz**2
        lower_band[:-1] += + A[+1:] / dz**2
        main_band[:] += -2 * A / dz**2

        ## Lower boundary condition (u=u_SSL)

        C[:nb] = u0
        main_band[:nb] = 1.
        upper_band[+1:nb+1] = 0.
        lower_band[+0:nb-1] = 0.

        if uf is None:
            ## Upper boundary condition (zero-gradient)

            # C[-nt:] = 0.
            # main_band[-nt:] = +1.
            # upper_band[-nt+1:] = +0.
            # lower_band[-nt-1:-1] = -2.
            # grad_band[-nt-2:-2] = +1.

            C[-nt:] = 0.
            upper_band[-nt+1:] = +0.
            main_band[-nt:] = +1.
            lower_band[-nt-1:-1] = -1.
            grad_band[-nt-2:-2] = +0.
        else:
            C[-nt:] = uf
            main_band[-nt:] = +1.
            upper_band[-nt+1:] = +0.
            lower_band[-nt-1:-1] = 0.
            grad_band[-nt-2:-2] = 0.

        # Assemble matrix

        M_banded = np.vstack([
            upper_band,
            main_band,
            lower_band,
            grad_band,
        ])

        # if not np.isfinite(M_banded).all():
        #     break

        np.savetxt("M_sparse.csv", M_banded, fmt="%10.2e", delimiter=",")
        np.savetxt("C_sparse.csv", C, fmt="%10.2e", delimiter=",")

        u = solve_banded((2, 1), M_banded, C)

        if bound_u0:
            u = np.maximum(u, u0)
        du = np.gradient(u, z)
        du = np.maximum(du, 0)

        if not np.isfinite(u).all():
            break

        if live_fig:

            if not plt.fignum_exists(fig.number):
                return u

            line.set_data(u, z)

            line_bis.set_data(u, z)
            ax_bis.relim()
            ax_bis.autoscale_view()

            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(1)
            print("*", iteration)

    plt.show(block=True)

    return u

z  = np.linspace(start, end, num=n, dtype=np.float64)
dz = (end - start) / (n - 1)

u = friedli(z, u_SSL, bound_u0=True)
print(f"{u = }")

# %%

def viscous_stress(u, du, e, de):
    return rho * nu * (de*u + e*du)

def turbulent_stress(u, du, Z_LS, e):
    phi = np.sqrt(Z_LS*u/nu)/A
    exponential_damping = 1 - np.exp(-phi)
    mixing_length = kappa * Z_LS * exponential_damping
    return + rho * e * mixing_length**2 * du**2

def dispersive_stress(u, du, e):
    mixing_length = dp * np.sqrt(lam * (1-e) /(1-eps_bulk))
    return + rho * e * mixing_length**2 * u/dp * du


def viscous_stress_rate(u, du, ddu, e, de, dde):
    return rho * nu * (dde*u + 2*de*du + e*ddu)

def turbulent_stress_rate(u, du, ddu, Z_LS, dZ_LS, e, de):
    phi = np.sqrt(Z_LS*u/nu) / A
    exponential_damping = 1 - np.exp(-phi)
    mixing_length = kappa * Z_LS * exponential_damping
    mixing_length_rate = kappa * (
        + dZ_LS * exponential_damping
        +  Z_LS * np.exp(-phi) * 1/2 * (dZ_LS*u + Z_LS*du) / np.sqrt(nu*A**2*Z_LS*u)
    )
    return + rho * (
        + de *                   mixing_length**2 *    du**2
        +  e * 2*mixing_length*mixing_length_rate *    du**2
        +  e *                   mixing_length**2 * 2*du*ddu
    )

def dispersive_stress_rate(u, du, ddu, e, de):
    mixing_length = dp * np.sqrt(lam * (1-e) /(1-eps_bulk))
    frac = np.zeros_like(e)
    mask = ~np.isclose(e, 1)
    frac[mask] = lam / ((1-eps_bulk) * (1-e[mask]))
    mixing_length_rate = - 1/2 * dp * np.sqrt(frac) * de
    return + rho / dp * (
        + de *                   mixing_length**2 *  u *  du
        +  e * 2*mixing_length*mixing_length_rate *  u *  du
        +  e *                   mixing_length**2 * du *  du
        +  e *                   mixing_length**2 *  u * ddu
    )


e = porosity(z)
du = np.gradient(u, z)
d2u = np.gradient(du, z)
de = porosity_rate(z)
d2e = porosity_accel(z)

e_s = porosity(z)
de_s = np.gradient(e_s, z)  # porosity_rate(z)
dde_s = np.gradient(de_s, z)  # porosity_accel(z)

Z_LS_s = integral_depth(z)
dZ_LS_s = np.gradient(Z_LS_s, z)  # integral_depth_rate(z)

tau_v_s = viscous_stress(u, du, e_s, de_s)
tau_t_s = turbulent_stress(u, du, Z_LS_s, e_s)
tau_d_s = dispersive_stress(u, du, e_s)
f = kozeny_carman(u, e)

fig, axes = plt.subplots(ncols=3, sharey=True, figsize=(10, 5))

line, = axes[0].plot(u, z, '-', ms=2, label=rf"$q={np.sum(u*(z[1]-z[0])):.2f}~\mathrm{{m^2/s}}$")
axes[0].set_ylabel("$z$")
axes[0].set_xlabel("$u$", c=line.get_color())
line, = axes[0].twiny().plot(e_s, z, '--', label=r"$\epsilon$", c="C1")
line.axes.set_xlabel(r"$\epsilon$", c=line.get_color())
axes[0].legend()

axes[1].plot(tau_v_s, z, '-', label=r"$\tau_v$")
axes[1].plot(tau_t_s, z, '-', label=r"$\tau_t$")
axes[1].plot(tau_d_s, z, '-', label=r"$\tau_d$")
G = np.cumsum(rho*g*i*e * (z[1]-z[0]))
axes[1].plot(-G, z, '--', alpha=0.5, label="$G$")

axes[1].set_xlabel(r"$\tau$")
axes[1].legend()

# axes[1].dataLim.x0 = - rho*g*i*(end - start)
# axes[1].dataLim.x1 = + rho*g*i*(end - start)

axes[2].fill_betweenx(z, -e*rho*g*i, fc="k", label=r"$\varepsilon\rho g i$", alpha=0.1)

axes[2].plot(np.gradient(tau_v_s, z), z, '-', label=r"$\tau_v'$", alpha=0.5, c="C0")
axes[2].plot(np.gradient(tau_t_s, z), z, '-', label=r"$\tau_t'$", alpha=0.5, c="C1")
axes[2].plot(np.gradient(tau_d_s, z), z, '-', label=r"$\tau_d'$", alpha=0.5, c="C2")

axes[2].plot(viscous_stress_rate(u, du, d2u, e_s, de_s, dde_s), z, '--', c="C0")
axes[2].plot(turbulent_stress_rate(u, du, d2u, Z_LS_s, dZ_LS_s, e_s, de_s), z, '--', c="C1")
axes[2].plot(dispersive_stress_rate(u, du, d2u, e_s, de_s), z, '--', c="C2")
axes[2].plot(e*f, z, ':', label=r"$\epsilon f$")

axes[2].plot(
    + viscous_stress_rate(u, du, d2u, e_s, de_s, dde_s)
    + turbulent_stress_rate(u, du, d2u, Z_LS_s, dZ_LS_s, e_s, de_s)
    + dispersive_stress_rate(u, du, d2u, e_s, de_s)
    + e*f
, z, 'k', zorder=0, lw=2, alpha=0.5, label=r"$\tau'+\epsilon f$")

axes[2].dataLim.x0 = -rho*g*i
axes[2].dataLim.x1 = +rho*g*i

axes[2].set_xlabel(r"contributions")
axes[2].legend()

plt.show()