# Two-dimensional  incompressible irrotational flow past a cylinder of radius a

using CairoMakie

a = 1.0 # radius of cylinder
U = 1.0 # free stream velocity

w(z, a, U) = U * (z + a^2 / z) # complex potential

xs = range(-2.5, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im, a, U)) for x in xs, y in ys]

contour(xs, ys, z, levels=range(-20, 20))