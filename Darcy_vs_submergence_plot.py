import numpy as np
from matplotlib import pyplot as plt

from fenics_attempt import VelocityProfileSolver

solver = VelocityProfileSolver()

solver.update({
    # Variables of the case
    "rho": 1000.0,
    "nu": 1.0e-6,
    "slope": 0.1e-2,
    "d_p": 8.0e-3,
    "eps_b": 0.38,
    "z0": -20.0e-3,
    "z1": 12.0e-3,
    "n_cells": 10_000,
    "lam": 0.1,
})

def ferguson(xi, a1=6.5, a2=2.5):
    C_f = 3 * xi / np.sqrt(1 + 0.18 * xi**1.67)  # Ferguson (2021)
    return 8 / C_f**2
    C_f = a1*a2*xi / np.sqrt(a1**2 + a2**2 * xi**(5/3))  # Eq. 1.83 (hydraulique)

def recking(xi, i=solver.slope):
    alpha = np.clip(4*xi**(-0.43), 1, 4) * np.clip(7*i**0.85*xi, 1, 2.6)
    return 8 / (6.25 + 5.75*np.log10(xi/alpha))**2  # Eq. 1.149 (hydraulique)
    C_f = 4.416 * xi**1.904 * (1 + (xi/1.283)**1.618)**(-1.083)
    return 8 / C_f**2

def gms(xi):
    return 8*9.81 / (23.1**2 * xi**(1/3))

def keulegan(xi):
    return (2.03*np.log10(12.2*xi/2.5))**-2  # k_s = 2.5 d_m, Eq. 1.73 (hydraulique)

num = 50

f = np.full(num, np.nan, dtype=np.float64)
# The datum sits at z = -d_p, so z1 = (xi - 1)*d_p and xi <= 1 puts the free
# surface at or below the grain crest: nothing to solve, and graded_mesh folds.
xi = np.logspace(np.log10(0.6), 2, num=num)

fig, (ax_u, ax_f) = plt.subplots(ncols=2)
l_f, = ax_f.plot([], [], '-o', mfc="none")
l_u, = ax_u.plot([], [])
# patch = ax_u.fill_betweenx([], [])
ax_f.plot(xi, gms(xi), label="Manning")
ax_f.plot(xi, ferguson(xi), label="Ferguson")
ax_f.plot(xi, recking(xi), label="Recking")
ax_f.plot(xi, keulegan(xi), label="Keulegan")
ax_f.loglog()
ax_f.legend()
ax_u.set_xlabel(r"$u/u_\ast$")
ax_u.set_ylabel(r"$\xi=h/d_p$")
ax_f.set_xlabel(r"$\xi=h/d_p$")
ax_f.set_ylabel(r"$f=8(u_\ast/\bar{u})^2$")
ax_f.yaxis.set_label_position("right")
ax_f.yaxis.tick_right()
fig.show()

z0 = 0.0 - solver.d_p / 2

for i, xii in enumerate(xi):

    solver.update({"z1": xii * solver.d_p + z0})
    zs, us = solver.run()

    # z0 = 0.
    depth = solver.z1 - z0
    # depth = solver.integral_depth(solver.z1)

    u_star = np.sqrt(depth * solver.g * solver.slope)
    m = zs > z0
    # `depth`, not `solver.depth`: the property is z1 - z(eps=0.8), a different
    # height from the z1 - z0 used for u_star and xi just above.  Mixing them
    # made u_mean too large by 2.1x at xi = 2, i.e. f 4.5x too small, with the
    # error worst at the low-submergence end.
    u_mean = np.trapezoid(us[m] * solver.porosity(zs[m]), zs[m]) / depth

    print(f"{u_mean / u_star = :g}  |  {depth / solver.d_p = :g}")
    f[i] = 8 * (u_star / u_mean)**2
    xi[i] = depth / solver.d_p

    if False:
        l_u.set_data(us / u_star, zs / solver.d_p)
    else:
        ax_u.plot(us / u_star, zs / solver.d_p)
    # patch.remove()
    # patch = ax_u.fill_betweenx(zs[m]/solver.d_p, us[m]/u_star, alpha=0.2)
    l_f.set_data(xi, f)
    for ax in [ax_u, ax_f]:
        ax.relim()
        ax.autoscale_view()
    fig.canvas.draw()
    fig.canvas.flush_events()


# ax_f.plot(xi, f**(0.8), '-o')
fig.canvas.draw()
fig.canvas.flush_events()
fig.tight_layout()
plt.show()

np.savetxt(
    "communications/hyporheic_flow/media/darcy_double_average.csv",
    np.column_stack((xi, f, gms(xi), recking(xi), keulegan(xi), ferguson(xi))),
    delimiter=",",
    header="xi,f,Manning,Recking,Keulegan,Ferguson",
    comments=""
)
