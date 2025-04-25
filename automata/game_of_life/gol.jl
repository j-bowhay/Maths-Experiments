using Pkg
Pkg.activate(temp=true)
Pkg.add("Plots")

##
using Plots


function neighbours(world::Matrix{Bool}, row::Int, col::Int, r::Int=1)
    m, n = size(world)
    rows = @. mod1(row + (-r:r), m)
    cols = @. mod1(col + (-r:r), n)
    return sum(i != row || j != col ? world[i, j] : 0 for i in rows, j in cols)
end


function neighbours(world::Matrix{Bool})
    m, n = size(world)
    return [neighbours(world, i, j) for i in 1:m, j in 1:n]
end


willsurvive(cell::Bool, k::Int) = k == 3 || k == 2 && cell


function evolve!(world::Matrix{Bool})
    n = neighbours(world)
    for i in eachindex(world)
        world[i] = willsurvive(world[i], n[i])
    end
end


function plot_game(world, frames, fps, filename)
    anim = @animate for i in 1:frame
        heatmap(world; axis = nothing, border = :none, cbar = false, ratio = :equal, yflip=true)
        evolve!(world)
    end
    gif(anim, "automata/game_of_life/$(filename).gif"; fps = fps)
end

## Glider

world = zeros(Bool, 20, 20)
world[1:3, 3] .= true
world[2, 1] = true
world[2, 3] = true
world[3, 2] = true

plot_game(world, 200, 30, "glider")

## Random
world = rand(100, 100) .< 0.5
plot_game(world, 1000, 30, "random")

## Pulsar
world = zeros(Bool, 17, 17)
line = zeros(17)
line[5:7] .= 1
line[11:13] .= 1

for ind in [3,8,10,15]
    world[ind, :] .= line
    world[:, ind] .= line
end

plot_game(world, 100, 30, "pulsar")