from matplotlib import pyplot as plt

from momentum_ODE_solver import Parameters

# Variables of the case
p = Parameters(
    rho=1000.0,
    nu=1.0e-6,
    slope=0.1e-3,
    d_p=8.0e-3,
    eps_b=1e-5,
    z0=-5.0e-1,
    z1=30.0e-2,
    n_cells=10_000,
)
num = 20

plt.xlabel(r"$u$ [m/s]")
plt.ylabel(r"$\xi=h/d_p$ [-]")
plt.grid()

double = {"name": "double", "ls": "-", "kappa": 0.41, "lam": 0.15}
no_rey = {"name": r"$\kappa=0$", "ls": "-.", "kappa": 0., "lam": 0.15}
no_disp = {"name": r"$\lambda=0$", "ls": "--", "kappa": 0.41, "lam": 0.}
pure_laminar = {"name": r"$\kappa=\lambda=0$",
                "ls": ":", "kappa": 0, "lam": 0.}

for params in [double, no_rey, no_disp, pure_laminar]:

    for k, v in params.items():
        if k not in ["name", "ls"]:
            setattr(p, k, v)
    try:
        r = p.solve()
    except RuntimeError:
        continue

    plt.plot(
        r.u, r.z, label=rf"$\lambda={r.lam}$, $\kappa={r.kappa}$", ls=params["ls"])

plt.legend()
plt.show()
