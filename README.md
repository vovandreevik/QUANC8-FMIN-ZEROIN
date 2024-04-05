# QUANC8-FMIN-ZEROIN

The program calculates the luminosity of a blackbody and explores the sustainability of solutions in the temperature range from $T = 1000°K$ to $T = 9000°K$ in increments of $T = 1000°K$

Using the FMIN program, the program calculates the value of $z^*$ by minimizing the function $f(z)=e^z (2z^2-4)+ (2z^2-1)^2+ e^2z- 3z^4$ in the interval $[-2,-1]$.

Multiplying $z^*$ by $(-3.039830×10^{-5})$ and getting $λ_1$ (the lower limit of integration).

Using the ZEROIN program, the program calculates the $y^*$ root of the equation: $2\sqrt{x}=cos \displaystyle\frac{πx}{2}$

Multiplying $y^*$ by $(31.66675×10^{-5})$ and getting $λ_2$ (the upper limit of integration).

Using the QUANC8 program, the luminosity is calculated (as a percentage) using the formula $EFF= \displaystyle\frac{64.77}{T^4}\int_{λ_1}^{λ_2} \displaystyle\frac{\mathrm{d}x}{x^5(e^{\displaystyle\frac{1.432}{Tx}} - 1)}$ in the temperature range from $T = 1000°K$ to $T = 9000°K$ in increments of $1000°$K.

The impact of errors in changing $λ_1$ and $λ_2$ on accuracy has been evaluated.

_`quanc8.h` and `Forsythe.h` are library programs_
