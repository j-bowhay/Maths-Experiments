# Two-dimensional  incompressible irrotational flow with a line vortex in a semi
# infinite channel

using CairoMakie

const Γ = 1
const d = 1.5
const a= 1
d_hat = sinh(π*d/2a)

w(z) = (im*Γ/2π)*(-log(sinh(π*z/2a) - d_hat) + log(sinh(π*z/2a) + d_hat))

xs = range(-1, 2.5, 200)
ys = range(-2.5, 2.5, 200)

z = [imag(w(x + y*im)) for x in xs, y in ys]

contour(xs, ys, z, levels=50)