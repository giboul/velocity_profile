import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import newton

from momentum_ODE_solver import Parameters, porosity

# Variables of the case
p = Parameters(
    rho=1000.0,
    nu=1.0e-6,
    slope=0.1e-2,
    d_p=8.0e-3,
    eps_b=0.38,
    z0=-5.0e-1,
    z1=12.0e-3,
    n_cells=10_000,
    lam=0.,
)
num = 10

xi = np.logspace(np.log10(0.6), 6, num=num)
u_star = np.full(num, np.nan, dtype=np.float64)
q_star = np.full(num, np.nan, dtype=np.float64)


def gms(q_star, k_manning=7):
    return k_manning**(3/5) * q_star**(2/5)


def _keulegan_objective_function(u_star, q_star, kappa=0.41, B=5.):
    return u_star - q_star**(1/3) * ((np.log(q_star/(2*u_star)) + 1) / kappa + B)**(2/3)


def keulegan(q_star, **fsolve_kwargs):
    return newton(_keulegan_objective_function, gms(q_star), args=[q_star], **fsolve_kwargs)


fig, (ax_u, ax_f) = plt.subplots(ncols=2)
l_f, = ax_f.plot([], [], '-o', mfc="none")
l_u, = ax_u.plot([], [])
ax_f.plot(_q_star := np.logspace(-3, 9), gms(_q_star), label="Manning")
ax_f.plot(_q_star := np.logspace(-3, 9), keulegan(_q_star), label="Keulegan")
ax_f.loglog()
ax_f.legend()
ax_u.set_xlabel(r"$u$ [m/s]")
ax_u.set_ylabel(r"$\xi=h/d_p$ [-]")
ax_f.set_xlabel(r"$q_\ast$ [-]")
ax_f.set_ylabel(r"$\bar{u}_\ast$")
ax_f.yaxis.set_label_position("right")
ax_f.yaxis.tick_right()
ax_f.grid()
ax_u.grid()
fig.show()

z0 = 0.0 - p.d_p / 4

for i, xii in enumerate(xi):

    if not plt.fignum_exists(fig.number):
        break

    p.z1 = xii * p.d_p + z0
    try:
        r = p.solve()
    except RuntimeError:
        break

    q = np.trapezoid(r.u * porosity(r.z, p), r.z)
    u_mean = (r.u * porosity(r.z, p))[p.z1 > p.z0].mean()

    u_star[i] = u_mean / np.sqrt(p.g * p.slope * p.d_p)
    q_star[i] = q / np.sqrt(p.g * p.slope * p.d_p**3)
    print(
        f"{u_star[i] = :.3e}   {q_star[i] = :.3e} | e = {r.total_momentum_error():.3g} %")

    ax_u.plot(r.u, r.z)
    l_f.set_data(q_star, u_star)
    for ax in [ax_u, ax_f]:
        ax.relim()
        ax.autoscale_view()
    fig.canvas.draw()
    fig.canvas.flush_events()


fig.canvas.draw()
fig.canvas.flush_events()
fig.tight_layout()
plt.show()
