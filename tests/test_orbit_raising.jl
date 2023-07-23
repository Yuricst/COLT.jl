"""
Test orbit raising problem
"""

using Plots 
using JuMP
using Memoize

include("../src/COLT.jl")

# construct problem
prob = COLT.CollocationProblem(COLT.get_orbit_raising_model(N=40)...)

# solve problem
COLT.solve!(prob)

# get model out to plot etc.
model = prob.model

@show objective_value(model);
@show termination_status(model);


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

initial_orbit = get_circle(1.0);
final_orbit   = get_circle(objective_value(model));

ts_controls = JuMP.value.(model[:t])
ts_states = JuMP.value.(model[:t])

rs = JuMP.value.(model[:r])
θs = JuMP.value.(model[:θ])
vrs = JuMP.value.(model[:vr])
vθs = JuMP.value.(model[:vθ])
u1s = JuMP.value.(model[:u1])
u2s = JuMP.value.(model[:u2])

xys = Vector[]
for (_r, _θ) in zip(rs,θs)
    push!(xys, [_r*cos(_θ), _r*sin(_θ)])
end
xys = hcat(xys...)

psh = plot(frame_style=:box, gridalpha=0.5, legend=:topleft)
plot!(psh, ts_states, rs, label="r", color=:red)
plot!(psh, ts_states, θs, label="θ", color=:blue)
plot!(psh, ts_states, vrs, label="vr", color=:purple)
plot!(psh, ts_states, vθs, label="vθ", color=:green)

pcon = plot(frame_style=:box, gridalpha=0.5, legend=:bottomleft,
    ylim=[-1.05, 1.05])
plot!(pcon, ts_controls, u1s, label="u1", color=:red, marker=:circle)
plot!(pcon, ts_controls, u2s, label="u2", color=:blue, marker=:circle)

ptraj = plot(frame_style=:box, aspect_ratio=:equal, legend=:bottomright)
plot!(ptraj, xys[1,:], xys[2,:], marker=:circle, label="Transfer", color=:lime)
plot!(ptraj, initial_orbit[1,:], initial_orbit[2,:], color=:blue, label="Initial")
plot!(ptraj, final_orbit[1,:], final_orbit[2,:], color=:red, label="Final")

l = @layout [[a;b] c]
pcomb = plot(psh, pcon, ptraj; size=(900,500), layout=l)
savefig(pcomb, joinpath(@__DIR__, "orbit_raising.png"))
display(pcomb)