"""
Test orbit raising problem
"""

using GLMakie 
using JuMP

include(joinpath(@__DIR__, "../src/COLT.jl"))

# construct problem
prob = COLT.HermiteSimpsonProblem(COLT.get_orbit_raising_model(N=40)...)

# solve problem
COLT.solve!(prob)

@show objective_value(prob.model);
@show termination_status(prob.model);


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
final_orbit   = get_circle(objective_value(prob.model));

ts_controls = JuMP.value.(prob.model[:t])
ts_states = JuMP.value.(prob.model[:t])

rs = JuMP.value.(prob.model[:r])
θs = JuMP.value.(prob.model[:θ])
vrs = JuMP.value.(prob.model[:vr])
vθs = JuMP.value.(prob.model[:vθ])
u1s = JuMP.value.(prob.model[:u1])
u2s = JuMP.value.(prob.model[:u2])

xys = Vector[]
for (_r, _θ) in zip(rs,θs)
    push!(xys, [_r*cos(_θ), _r*sin(_θ)])
end
xys = hcat(xys...)

fig = Figure(resolution = (800, 500))
ax1 = Axis(fig[1:2, 1], xlabel = "x", ylabel = "y", aspect = 1)
scatterlines!(ax1, xys[1,:], xys[2,:], label="Transfer", color=:grey22, marker=:circle, markersize=7.0)
lines!(ax1, initial_orbit[1,:], initial_orbit[2,:], color=:blue, label="Initial")
lines!(ax1, final_orbit[1,:], final_orbit[2,:], color=:red, label="Final")
axislegend(ax1, position=:cc)

ax2 = Axis(fig[1, 2], xlabel = "t", ylabel = "State")
scatterlines!(ax2, ts_states, rs, label="r", color=:red, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts_states, θs, label="θ", color=:blue, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts_states, vrs, label="vr", color=:purple, marker=:circle, markersize=7.0)
scatterlines!(ax2, ts_states, vθs, label="vθ", color=:lime, marker=:circle, markersize=7.0)
axislegend(ax2, position=:lt)

ax3 = Axis(fig[2, 2], xlabel = "t", ylabel = "Control")
scatterlines!(ax3, ts_controls, u1s, label="u1", color=:red, marker=:circle, markersize=7.0)
scatterlines!(ax3, ts_controls, u2s, label="u2", color=:blue, marker=:circle, markersize=7.0)
axislegend(ax3, position=:lb)
save(joinpath(@__DIR__, "orbit_raising.png"), fig)
fig

# psh = plot(frame_style=:box, gridalpha=0.5, legend=:topleft)
# plot!(psh, ts_states, rs, label="r", color=:red)
# plot!(psh, ts_states, θs, label="θ", color=:blue)
# plot!(psh, ts_states, vrs, label="vr", color=:purple)
# plot!(psh, ts_states, vθs, label="vθ", color=:green)

# pcon = plot(frame_style=:box, gridalpha=0.5, legend=:bottomleft,
#     ylim=[-1.05, 1.05])
# plot!(pcon, ts_controls, u1s, label="u1", color=:red, marker=:circle)
# plot!(pcon, ts_controls, u2s, label="u2", color=:blue, marker=:circle)

# ptraj = plot(frame_style=:box, aspect_ratio=:equal, legend=:bottomright)
# plot!(ptraj, xys[1,:], xys[2,:], marker=:circle, label="Transfer", color=:lime)
# plot!(ptraj, initial_orbit[1,:], initial_orbit[2,:], color=:blue, label="Initial")
# plot!(ptraj, final_orbit[1,:], final_orbit[2,:], color=:red, label="Final")

# l = @layout [[a;b] c]
# pcomb = plot(psh, pcon, ptraj; size=(900,500), layout=l)
# savefig(pcomb, joinpath(@__DIR__, "orbit_raising.png"))
# display(pcomb)