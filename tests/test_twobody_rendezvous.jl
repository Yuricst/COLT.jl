"""
Test twobody rendez-vous problem
"""

using Plots 
using JuMP
using Memoize

include("../src/COLT.jl")

# construct problem
mu = 1.0
c1 = 1e-3
c2 = 1e-3
state0 = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0]
target_state = [1.2, 0.0, 0.0, 0.0, sqrt(1/1.2), 0.01]
tf_bounds = [0.3, 3.0]

prob = COLT.CollocationProblem(COLT.get_twobody_rendezvous_model(
    state0,
    target_state,
    mu,
    c1,
    c2,
    tf_bounds;
)...)

# solve problem
COLT.solve!(prob)

# get model out to plot etc.
model = prob.model

tf = JuMP.value.(model[:tf])
xs = JuMP.value.(model[:x])
ys = JuMP.value.(model[:y])
zs = JuMP.value.(model[:z])

println("tf = ", tf)

ptraj = plot(frame_style=:box, aspect_ratio=:equal, legend=:bottomright)
plot!(ptraj, xs, ys)
display(ptraj)