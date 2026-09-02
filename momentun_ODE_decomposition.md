# Governing equation

$$0 = \epsilon\varrho g i + \epsilon\bar f + \left(\tau_v + \tau_t + \tau_d\right)' - \mu\epsilon'u'$$

## Viscous stress rate

$$\tau_v = \mu\left(\epsilon u\right)' \Rightarrow \tau_v'=\mu\epsilon''u+2\mu\epsilon'u'+\mu\epsilon u''$$

## Turbulent stress rate

$$\tau_t = -\varrho\epsilon\kappa^2 Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2$$

$$\begin{align*}
\Rightarrow -\frac{\tau_t'}{\varrho\kappa^2} = & + \epsilon' Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2\\
    & + 2 \epsilon Z_{LS} Z_{LS}' \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2\\
    & + \epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}'u + Z_{LS}u'}{\sqrt{\nu A_\ast^2 Z_{LS}u}} u'^2\\
    & + 2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u' u''\\
\end{align*}$$

## Dispersive stress rate

$$\tau_d = - \frac{\varrho\lambda d_p}{1-\epsilon_b}[1-\epsilon]\epsilon uu'$$

$$\begin{align*}
\Rightarrow -\frac{1-\epsilon_b}{\varrho \lambda d_p}\tau_d' = & + (1-2\epsilon)\epsilon' u u'\\
    & + (1-\epsilon)\epsilon u'^2\\
    & + (1-\epsilon)\epsilon u u''\\
\end{align*}$$

# Expanded governing equation

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & + \mu(\epsilon u)'' - \mu u'\epsilon'\\
    & -\varrho\kappa^2\left(\epsilon Z_{LS}^2\left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2u'^2\right)'\\
    & - \frac{\varrho d_p\lambda}{1-\epsilon_b}\left([1-\epsilon]\epsilon uu'\right)'\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \frac{B_E(1-\epsilon)\epsilon}{d_p}\varrho u^2
\end{align*}$$

$$\begin{align*}
0 = & + \epsilon\varrho gi\\
    & + \mu\epsilon u'' + \mu u'\epsilon' + \mu\epsilon'' u\\
    & - \varrho \kappa^2 \left[ \epsilon' Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + 2 \epsilon Z_{LS} Z_{LS}' \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + \epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}'u + Z_{LS}u'}{\sqrt{\nu A_\ast^2 Z_{LS}u}} u'^2 + 2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u' u'' \right]\\
    & -\frac{\varrho \lambda d_p}{1-\epsilon_b} \left[ (1-2\epsilon)\epsilon' u u' + (1-\epsilon)\epsilon u'^2 + (1-\epsilon)\epsilon u u'' \right]\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \frac{B_E(1-\epsilon)\epsilon}{d_p}\varrho u^2
\end{align*}$$

## Regrouping $u$, $u'$, $u'^2$, $u'^3$ and $u''$

$$\begin{align*}
0 = & + \epsilon\varrho gi + \mu\epsilon'' u - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \frac{B_E(1-\epsilon)\epsilon}{d_p}\varrho u^2\\
    & + \left[\mu \epsilon' -\frac{\varrho \lambda d_p}{1-\epsilon_b}(1-2\epsilon)\epsilon' u \right] u'\\
    & - \left\{\varrho \kappa^2\left[\epsilon' Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 + 2 \epsilon Z_{LS} Z_{LS}' \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 + \epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}'u}{\sqrt{\nu A_\ast^2 Z_{LS}u}} \right] + \frac{\varrho \lambda d_p}{1-\epsilon_b} (1-\epsilon)\epsilon\right\} u'^2\\
    & - \varrho \kappa^2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}}{\sqrt{\nu A_\ast^2 Z_{LS}u}} u'^3 \\
    & + \left[\mu\epsilon - \frac{\varrho \lambda d_p}{1-\epsilon_b}(1-\epsilon)\epsilon u \right] u''\\
    & - 2 \varrho \kappa^2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'  u''\\
\end{align*}$$



# Expressing the second derivative

$$\begin{align*}
 - \left[ \mu\epsilon - 2 \varrho \kappa^2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'  -\frac{\varrho \lambda d_p}{1-\epsilon_b}(1-\epsilon)\epsilon u \right] u'' = &
 + \epsilon\varrho gi\\
    & + \mu u'\epsilon' + \mu\epsilon'' u\\
    & - \varrho \kappa^2 \left[ \epsilon' Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + 2 \epsilon Z_{LS} Z_{LS}' \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + \epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}'u + Z_{LS}u'}{\sqrt{\nu A_\ast^2 Z_{LS}u}} u'^2 \right]\\
    & -\frac{\varrho \lambda d_p}{1-\epsilon_b} \left[ (1-2\epsilon)\epsilon' u u'+ (1-\epsilon)\epsilon u'^2\right]\\
    & - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \frac{B_E(1-\epsilon)\epsilon}{d_p}\varrho u^2
\end{align*}$$


$$u'' = -
\frac{
    \epsilon\varrho gi + \mu u'\epsilon' + \mu\epsilon'' u - \varrho \kappa^2 \left[ \epsilon' Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + 2 \epsilon Z_{LS} Z_{LS}' \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'^2 + \epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right) e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}} \frac{Z_{LS}'u + Z_{LS}u'}{    \sqrt{\nu A_\ast^2 Z_{LS}u}} u'^2 \right] -\frac{\varrho \lambda d_p}{1-\epsilon_b} \left[ (1-2\epsilon)\epsilon' u u'+ (1-\epsilon)\epsilon u'^2\right] - \frac{A_E(1-\epsilon)^2}{\epsilon d_p^2}\mu u - \frac{B_E(1-\epsilon)\epsilon}{d_p}\varrho u^2
}{
    \mu\epsilon - 2 \varrho \kappa^2\epsilon Z_{LS}^2 \left(1-e^{-\sqrt{\frac{Z_{LS}u}{\nu A_\ast^2}}}\right)^2 u'  -\frac{\varrho \lambda d_p}{1-\epsilon_b}(1-\epsilon)\epsilon u
}$$