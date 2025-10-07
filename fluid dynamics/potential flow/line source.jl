# Two-dimensional  incompressible irrotational flow with a line source

using CairoMakie

const Q = 1

w(z) = (Q/2π)*log(z) # complex potential

xs = range(-2.5, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im)) for x in xs, y in ys]

contour(xs, ys, z, levels=10)