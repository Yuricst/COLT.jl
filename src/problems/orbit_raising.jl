"""
Ready-to-go model for orbit raising problem
"""


"""
Dynamics for orbit raising problem
"""
function dynamics_orbitraise(t, x, u, p)
    r, θ, vr, vθ = x  # unpack states
    a = 0.1405/(1 - 0.0749t)
    dx = [
        vr,
        vθ/r,
        vθ^2/r - 1/r^2 + a*u[1],
        -vθ*vr/r + a*u[2]
    ]
    return dx
end



"""
Get JuMP model and defect function for orbit raising problem
"""
function get_orbit_raising_model(;
    N::Int=30,
    solver=Ipopt.Optimizer,
    collocation_type="hermite_simpson"
)
    # problem parameters
    nx, nu = 4, 2         # number of states & controls
    tf_max = 3.32         # fixed final time
    N_nodes = N + 1       # number of nodes
    N_controls = N + 1    # number of controls

    # initialize model
    model = Model(solver; add_bridges = false)
    #set_optimizer_attribute(model, "algorithm", :LD_MMA)

    # decision variables
    @variables(model, begin
        # full duration of trajectory
        tf
        # times --- these are place holders, they will be fixed for now
        t[1:N_controls]
        # states
        0.5 ≤ r[1:N_nodes] ≤ 2.0
        0.0 ≤ θ[1:N_nodes] ≤ π
        0.0 ≤ vr[1:N_nodes] ≤ 2.0
        0.0 ≤ vθ[1:N_nodes] ≤ 2.0
        # controls
        -1.0 ≤ u1[1:N_controls] ≤ 1.0
        -1.0 ≤ u2[1:N_controls] ≤ 1.0
    end);

    # fix time
    fix(tf, tf_max; force = true)
    for i = 1:N_controls
        fix(t[i], (i-1)/(N_controls-1); force = true)
    end

    # boundary constraints (initial)
    fix(r[1],  1.0; force = true)
    fix(θ[1],  0.0; force = true)
    fix(vr[1], 0.0; force = true)
    fix(vθ[1], 1.0; force = true)

    # boundary constraints (final)
    fix(vr[end], 0.0; force = true)
    @constraint(model, sqrt(1/r[end]) - vθ[end] == 0.0)

    # path constraints
    for i = 1:N_controls
        @constraint(model, u1[i]^2 + u2[i]^2 ≤ 1.0)
    end

    # objective is linear, to maximize final radius
    @objective(model, Max, r[end]);

    # list out states and controls variables
    state_variables = [r,θ,vr,vθ]
    control_variables = [u1,u2]
    
    return model, dynamics_orbitraise, state_variables, control_variables
end