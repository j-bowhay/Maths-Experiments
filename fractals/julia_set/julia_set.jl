using Plots

L = 1500
K = 1000
x = LinRange(-1.5, 1.5, L)
y = LinRange(-1, 1, K)
c = -0.4 + 0.61im
R = 2
N = 1000

function juliaset(z, c, R, N)
    for n = 0:N
        if abs(z) > R^2 - R
            return n / N
        end
        z = z^2 + c
    end
    return 0
end

## Static version

A = @. juliaset(x' + y * im, c, R, N)

heatmap(A; c = :viridis, clims = (0, 0.15), cbar = :none, axis = :none, ticks = :none)

## Animation

cs = 0.7885 .* exp.(range(π / 2, 3π / 2; length = 50) .* im)
anim = @animate for c in cs
    A = juliaset.(x' .+ y .* im, c, R, N)
    heatmap(
        A;
        c = :viridis,
        clims = (0, 0.15),
        cbar = :none,
        axis = :none,
        ticks = :none,
        size = (800, 600),
    )
end
gif(anim, "fractals/julia_set/juliaset.gif", fps = 20)
