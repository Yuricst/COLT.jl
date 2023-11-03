"""
Test twobody rendez-vous problem
"""

using GLMakie 
using JuMP
using LinearAlgebra

include(joinpath(@__DIR__, "../../AstrodynamicsBase.jl/src/AstrodynamicsBase.jl"))
include(joinpath(@__DIR__, "../src/COLT.jl"))

# construct problem
mu = 1.0
c1 = 1e-3
c2 = 1e-3
state0_kep = [1.0, 0.0, 0.0, 0.0, 0.0, 0.3]
statef_kep = [1.5, 0.1, 0.2, 0.2, 0.1, 4.0]

state0 = vcat(AstrodynamicsBase.kep2cart(state0_kep, 1.0), [1.0])
target_state = AstrodynamicsBase.kep2cart(statef_kep, 1.0)
tf_bounds = [3.0, 8.0]

prob = COLT.HermiteSimpsonProblem(COLT.get_twobody_rendezvous_model(
    state0,
    target_state,
    mu,
    c1,
    c2,
    tf_bounds;
    N=50,
)...)

# solve problem
COLT.solve!(prob)

# get model out to plot etc.
model = prob.model

tf = JuMP.value.(model[:tf])
ts = JuMP.value.(model[:t])
xs = JuMP.value.(model[:x])
ys = JuMP.value.(model[:y])
zs = JuMP.value.(model[:z])
vxs = JuMP.value.(model[:vx])
vys = JuMP.value.(model[:vy])
vzs = JuMP.value.(model[:vz])
u1s = JuMP.value.(model[:u1])
u2s = JuMP.value.(model[:u2])
u3s = JuMP.value.(model[:u3])


# plotting
function get_circle(radius=1.0, center=[0,0], steps=100)
    xys = zeros(2,steps)
    thetas = LinRange(0,2π,steps)
    for i = 1:steps
        xys[1,i] = radius*cos(thetas[i]) + center[1]
        xys[2,i] = radius*sin(thetas[i]) + center[2]
    end
    return xys
end

initial_orbit = AstrodynamicsBase.get_orbit(state0_kep, 1.0, "keplerian")
final_orbit   = AstrodynamicsBase.get_orbit(statef_kep, 1.0, "keplerian")

fig = Figure(resolution = (800, 500))
ax1 = Axis(fig[1:2, 1], xlabel = "x", ylabel = "y", aspect = 1)
scatterlines!(ax1, xs, ys, label="Transfer", color=:grey22, marker=:circle, markersize=7.0)
lines!(ax1, initial_orbit[1,:], initial_orbit[2,:], color=:blue, label="Initial")
lines!(ax1, final_orbit[1,:], final_orbit[2,:], color=:red, label="Final")
axislegend(ax1, position=:cc)

ax2 = Axis(fig[1, 2], xlabel = "t", ylabel = "State")
scatterlines!(ax2, ts, xs, label="x", color=:red, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts, ys, label="y", color=:blue, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts, zs, label="z", color=:purple, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts, vxs, label="vx", color=:orange, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts, vys, label="vy", color=:green, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts, vzs, label="vz", color=:black, marker=:circle, markersize=7.0)
axislegend(ax2, position=:lt)

ax3 = Axis(fig[2, 2], xlabel = "t", ylabel = "Control")
scatterlines!(ax3, ts, u1s, label="u1", color=:red, marker=:circle, markersize=7.0)
scatterlines!(ax3, ts, u2s, label="u2", color=:blue, marker=:circle, markersize=7.0)
scatterlines!(ax3, ts, u3s, label="u3", color=:purple, marker=:circle, markersize=7.0)
axislegend(ax3, position=:lb)

save(joinpath(@__DIR__, "twobody_rendezvous.png"), fig)
fig