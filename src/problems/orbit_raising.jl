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
    model = Model(solver)
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
    @NLconstraint(model, sqrt(1/r[end]) - vθ[end] == 0.0)

    # path constraints
    for i = 1:N_controls
        @NLconstraint(model, u1[i]^2 + u2[i]^2 ≤ 1.0)
    end

    # objective is linear, to maximize final radius
    @objective(model, Max, r[end]);

    # list out states and controls variables
    state_variables = [r,θ,vr,vθ]
    control_variables = [u1,u2]

    if collocation_type == "hermite_simpson"
        hs_defect = get_hs_defect_function(
            dynamics_orbitraise,
            nx,
            nu,
            nothing,
        )
    else
        error("collocation_type $collocation_type not recognized")
    end

    # make function aliases with single scalar outputs
    hs_defect_1(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[1]
    hs_defect_2(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[2]
    hs_defect_3(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[3]
    hs_defect_4(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[4]

    # register & add nonlinear constraints for dynamics
    register(model, :hs_defect_1, 2*(1+nx+nu)+1, hs_defect_1; autodiff = true)
    register(model, :hs_defect_2, 2*(1+nx+nu)+1, hs_defect_2; autodiff = true)
    register(model, :hs_defect_3, 2*(1+nx+nu)+1, hs_defect_3; autodiff = true)
    register(model, :hs_defect_4, 2*(1+nx+nu)+1, hs_defect_4; autodiff = true)

    for i = 1:N
        ix1, ix2 = i, i+1
        iu1, iu2 = i, i+1
        txu1 = [t[ix1], r[ix1], θ[ix1], vr[ix1], vθ[ix1], u1[iu1], u2[iu1]]
        txu2 = [t[ix2], r[ix2], θ[ix2], vr[ix2], vθ[ix2], u1[iu2], u2[iu2]]
        @NLconstraint(model, hs_defect_1(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_2(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_3(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_4(txu1..., txu2..., tf) == 0.0)
    end

    # append objective
    @objective(model, Max, r[end]);
    return model, hs_defect
end