# Two-dimensional incompressible irrotational uniform flow past a flat plate

using CairoMakie

α = π/6
a = 1
U = 1.0 # free stream velocity
Γ = -4π*a*U*sin(α)

w(z) = U*z*cos(α) - sign(real(z))*im*U*sqrt(z^2 - 4a^2)*sin(α) - (im*Γ/(2π))*log((z + sign(real(z))*sqrt(z^2 - 4a^2))/2)

xs = range(-4, 4, 200)
ys = range(-4, 4, 200)

z = [imag(w(x + y*im)) for x in xs, y in ys]

fig, ax = contour(xs, ys, z, levels=40)
lines!(ax, [-2a, 2a], [0, 0], color=:black,)
fig