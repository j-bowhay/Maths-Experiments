# Two-dimensional incompressible irrotational uniform flow

using CairoMakie

α = π/4
U = 1.0 # free stream velocity

w(z, α, U) = U * exp(-α*im)*z # complex potential

xs = range(-2.5, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im, α, U)) for x in xs, y in ys]

contour(xs, ys, z, levels=range(-20, 20))