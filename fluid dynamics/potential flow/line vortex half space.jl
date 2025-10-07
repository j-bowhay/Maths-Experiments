# Two-dimensional  incompressible irrotational flow with a line vortex in half space

using CairoMakie

const Γ = 1
const d = 1.5

w(z) = (-im*Γ/2π)*log(z - d) + (im*Γ/2π)*log(z + d) # complex potential

xs = range(-1, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im)) for x in xs, y in ys]

contour(xs, ys, z, levels=30)