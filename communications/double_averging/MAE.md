---
title: Matched Asymptotic expansions for a uniform flow over a steep rough bed
author:
    - Axel Giboulot
date: \today
---
# Introduction

While multiple flow resistance laws and velocity profiles have been proposed, only a few account for hyporheic flow [@Lamb2007].

# Uniform, steady-state conditions

$$0 = \epsilon\varrho gi - \epsilon\bar f + \frac{\partial}{\partial z}\left(\tau_v + \tau_t + \tau_d\right) - \mu\frac{\partial u}{\partial z}\frac{\partial \epsilon}{\partial z}.$$

Developpig this exression yields

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & + \epsilon\bar f\\
    & + \partial_z(\mu\partial_z(\epsilon u)) - \mu\partial_z u\partial_z\epsilon\\
    & -\varrho\partial_z(\epsilon u'w')\\
    & -\varrho\partial_z(\epsilon \tilde u\tilde w)\\
\end{align*}$$

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \varrho\epsilon\frac{B_E(1-\epsilon)}{d_p}u^2\\
    & + \mu\left[\partial_z^2u\epsilon + \mu\partial_z u\partial_z\epsilon + \partial_z^2\epsilon u\right]\\
    & -\varrho\partial_z\left(\epsilon\ell_t^2(\partial_z u)^2\right)\\
    & -\varrho\partial_z\left(\epsilon\ell_d^2\frac{u}{d_p}\partial_z u\right)\\
\end{align*}$$

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \varrho\epsilon\frac{B_E(1-\epsilon)}{d_p}u^2\\
    & + \mu\left[\partial_z^2u\epsilon + \mu\partial_z u\partial_z\epsilon + \partial_z^2\epsilon u\right]\\
    & -\varrho\kappa^2\partial_z\left[\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2(\partial_z u)^2\right]\\
    & -\varrho d_p\lambda\partial_z\left(\epsilon\frac{1-\epsilon}{1-\epsilon_b} u\partial_z u\right)\\
\end{align*}$$

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \varrho\epsilon\frac{B_E(1-\epsilon)}{d_p}u^2\\
    & + \mu\left[\epsilon\partial_z^2u + \partial_z u\partial_z\epsilon + u\partial_z^2\epsilon\right]\\
    & -\varrho\kappa^2\left[\partial_z\left(\epsilon Z_{LS}^2\right) \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2(\partial_z u)^2 + 2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\sqrt{\frac{Z_{LS}}{\nu A_\ast^2u}} (\partial_z u)^2 + 2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2\partial_z u\partial_z^2u\right]\\
    & -\varrho d_p\lambda\left(\partial_z\left(\epsilon\frac{1-\epsilon}{1-\epsilon_b}\right) u\partial_z u + \epsilon\frac{1-\epsilon}{1-\epsilon_b} \left(\partial_z u\right)^2 + \epsilon\frac{1-\epsilon}{1-\epsilon_b} u\partial_z^2 u\right)\\
\end{align*}$$