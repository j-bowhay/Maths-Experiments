# Two-dimensional  incompressible irrotational stagnation point flow

using CairoMakie

w(z) = z^2/2 # complex potential

xs = range(-2.5, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im)) for x in xs, y in ys]

contour(xs, ys, z, levels=range(-20, 20))