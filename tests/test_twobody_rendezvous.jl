"""
Test twobody rendez-vous problem
"""

using Plots 
using JuMP
using LinearAlgebra
using Memoize
import AstrodynamicsBase

include("../src/COLT.jl")

# construct problem
mu = 1.0
c1 = 1e-3
c2 = 1e-3
state0_kep = [1.0, 0.0, 0.0, 0.0, 0.0, 0.3]
statef_kep = [1.5, 0.1, 0.2, 0.2, 0.1, 0.9]

state0 = vcat(AstrodynamicsBase.kep2cart(state0_kep, 1.0), [1.0])
target_state = AstrodynamicsBase.kep2cart(statef_kep, 1.0)
tf_bounds = [3.0, 8.0]

prob = COLT.CollocationProblem(COLT.get_twobody_rendezvous_model(
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

println("tf = ", tf)


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

ptraj = plot(frame_style=:box, aspect_ratio=:equal, legend=:bottomright)
plot!(ptraj, xs, ys, marker=:circle, label="Transfer", color=:lime)
plot!(ptraj, initial_orbit[1,:], initial_orbit[2,:], color=:blue, label="Initial")
plot!(ptraj, final_orbit[1,:], final_orbit[2,:], color=:red, label="Final")

psh = plot(frame_style=:box, gridalpha=0.5, legend=:topleft)
plot!(psh, ts, xs, label="x", marker=:circle)
plot!(psh, ts, ys, label="y", marker=:circle)
plot!(psh, ts, zs, label="z", marker=:circle)
plot!(psh, ts, vxs, label="vx", marker=:circle)
plot!(psh, ts, vys, label="vy", marker=:circle)
plot!(psh, ts, vzs, label="vz", marker=:circle)

pcon = plot(frame_style=:box, gridalpha=0.5, legend=:bottomleft,
    ylim=[-1.05, 1.05])
plot!(pcon, ts, u1s, label="u1", marker=:circle)
plot!(pcon, ts, u2s, label="u2", marker=:square)
plot!(pcon, ts, u3s, label="u3", marker=:diamond)

l = @layout [[a;b] c]
pcomb = plot(psh, pcon, ptraj; size=(900,500), layout=l)
display(pcomb)