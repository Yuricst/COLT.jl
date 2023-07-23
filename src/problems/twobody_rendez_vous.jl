"""
Ready-to-go model for two-body spacecraft rendez-vous
"""


"""
Dynamics for two-body spacecraft rendez-vous
"""
function dynamics_twobody(t, states, u, p)
    x, y, z, _, _, _, mass = states  # unpack states
    r3 = (x^2 + y^2 + z^2)^(3/2)    # radius
    mu, c1, c2 = p                   # unpack parameters
    # dx = [
    #     states[4],
    #     states[5],
    #     states[6],
    #     -mu/r3*x + c1/mass*u[1],
    #     -mu/r3*y + c1/mass*u[2],
    #     -mu/r3*z + c1/mass*u[3],
    #     -c2*sqrt(u[1]^2 + u[2]^2 + u[3]^2)
    # ]
    dx = [
        states[4],
        states[5],
        states[6],
        -mu/r3*x + c1/mass*u[1],
        -mu/r3*y + c1/mass*u[2],
        -mu/r3*z + c1/mass*u[2],
        -c2*(u[1]^2 + u[2]^2 + u[3]^2)^(0.5)
    ]
    return dx
end


"""
Get JuMP model and defect function for two-body spacecraft rendez-vous
"""
function get_twobody_rendezvous_model(
    state0,
    target_state,
    mu, c1, c2, 
    tf_bounds,
    mass_bounds=[0.3, 1.0],
    r_max::Real=3.0,
    v_max::Real=1.5,
    N::Int=40,
    solver=Ipopt.Optimizer,
    collocation_type="hermite_simpson",
    tf_free::Bool=true,
    max_mf::Bool=true,
)
    # problem parameters
    nx, nu = 7, 3         # number of states & controls
    N_nodes = N + 1       # number of nodes
    N_controls = N + 1    # number of controls

    # initialize model
    model = Model(solver)

    # decision variables
    @variables(model, begin
        # full duration of trajectory
        tf_bounds[1] ≤ tf ≤ tf_bounds[2]
        # times --- these are place holders, they will be fixed for now
        t[1:N_controls]
        # states
        -r_max ≤ x[1:N_nodes] ≤ r_max, (start = 1.0)
        -r_max ≤ y[1:N_nodes] ≤ r_max, (start = 0.1)
        -r_max ≤ z[1:N_nodes] ≤ r_max, (start = 0.1)
        -v_max ≤ vx[1:N_nodes] ≤ v_max, (start = 0.1)
        -v_max ≤ vy[1:N_nodes] ≤ v_max, (start = 1.1)
        -v_max ≤ vz[1:N_nodes] ≤ v_max, (start = 0.1)
        mass_bounds[1] ≤ mass[1:N_nodes] ≤ mass_bounds[2], (start = 1.0)
        # controls
        -1.0 ≤ u1[1:N_controls] ≤ 1.0
        -1.0 ≤ u2[1:N_controls] ≤ 1.0
        -1.0 ≤ u3[1:N_controls] ≤ 1.0
    end);

    # fix final time
    if tf_free == false
        fix(tf, tf_bounds[2]; force = true)
    end

    # fix node locations
    for i = 1:N_controls
        fix(t[i], (i-1)/(N_controls-1); force = true)   # fix the node locations 
    end

    # boundary constraints (initial)
    fix(x[1],  state0[1]; force = true)
    fix(y[1],  state0[2]; force = true)
    fix(z[1],  state0[3]; force = true)
    fix(vx[1], state0[4]; force = true)
    fix(vx[1], state0[5]; force = true)
    fix(vx[1], state0[6]; force = true)
    fix(mass[1], state0[7]; force = true)

    # boundary constraints (final)
    if tf_free == false
        # compute final targeted state
        fix(x[end],  target_state[1]; force = true)
        fix(y[end],  target_state[2]; force = true)
        fix(z[end],  target_state[3]; force = true)
        fix(vx[end], target_state[4]; force = true)
        fix(vy[end], target_state[5]; force = true)
        fix(vz[end], target_state[6]; force = true)
    else
        # add NL constraint based on rendez-vous
        @memoize get_target_state(_t) = keplerder_nostm(mu, target_state, 0.0, _t, 1.e-12, 10)
        get_target_state_1(t) = get_target_state(t)[1]
        get_target_state_2(t) = get_target_state(t)[2]
        get_target_state_3(t) = get_target_state(t)[3]
        get_target_state_4(t) = get_target_state(t)[4]
        get_target_state_5(t) = get_target_state(t)[5]
        get_target_state_6(t) = get_target_state(t)[6]
        register(model, :get_target_state_1, 1, get_target_state_1; autodiff = true)
        register(model, :get_target_state_2, 1, get_target_state_2; autodiff = true)
        register(model, :get_target_state_3, 1, get_target_state_3; autodiff = true)
        register(model, :get_target_state_4, 1, get_target_state_4; autodiff = true)
        register(model, :get_target_state_5, 1, get_target_state_5; autodiff = true)
        register(model, :get_target_state_6, 1, get_target_state_6; autodiff = true)
    end
    @NLconstraint(model, get_target_state_1(tf) - x[end] == 0.0)
    @NLconstraint(model, get_target_state_2(tf) - y[end] == 0.0)
    @NLconstraint(model, get_target_state_3(tf) - z[end] == 0.0)
    @NLconstraint(model, get_target_state_4(tf) - vx[end] == 0.0)
    @NLconstraint(model, get_target_state_5(tf) - vy[end] == 0.0)
    @NLconstraint(model, get_target_state_6(tf) - vz[end] == 0.0)

    # path constraints
    for i = 1:N_controls
        @NLconstraint(model, u1[i]^2 + u2[i]^2 + u3[i]^2 ≤ 1.0)
    end

    # objective is linear, to maximize final mass or minimize time
    if max_mf == true
        @objective(model, Max, mass[end]);
    else
        @objective(model, Min, tf);
    end

    if collocation_type == "hermite_simpson"
        hs_defect = get_hs_defect_function(
            dynamics_twobody,
            nx,
            nu,
            [mu, c1, c2],
        )
    else
        error("collocation_type $collocation_type not recognized")
    end

    # make function aliases with single scalar outputs
    hs_defect_1(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[1]
    hs_defect_2(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[2]
    hs_defect_3(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[3]
    hs_defect_4(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[4]
    hs_defect_5(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[5]
    hs_defect_6(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[6]
    hs_defect_7(txu0_txu1_uc...) = hs_defect(txu0_txu1_uc...)[7]

    # register & add nonlinear constraints for dynamics
    register(model, :hs_defect_1, 2*(1+nx+nu)+1, hs_defect_1; autodiff = true)
    register(model, :hs_defect_2, 2*(1+nx+nu)+1, hs_defect_2; autodiff = true)
    register(model, :hs_defect_3, 2*(1+nx+nu)+1, hs_defect_3; autodiff = true)
    register(model, :hs_defect_4, 2*(1+nx+nu)+1, hs_defect_4; autodiff = true)
    register(model, :hs_defect_5, 2*(1+nx+nu)+1, hs_defect_5; autodiff = true)
    register(model, :hs_defect_6, 2*(1+nx+nu)+1, hs_defect_6; autodiff = true)
    register(model, :hs_defect_7, 2*(1+nx+nu)+1, hs_defect_7; autodiff = true)

    for i = 1:N
        ix1, ix2 = i, i+1
        iu1, iu2 = i, i+1
        txu1 = [t[ix1], x[ix1], y[ix1], z[ix1],
                vx[ix1], vy[ix1], vz[ix1], mass[ix1],
                u1[iu1], u2[iu1], u3[iu1]]

        txu2 = [t[ix2], x[ix2], y[ix2], z[ix2],
                vx[ix2], vy[ix2], vz[ix2], mass[ix2],
                u1[iu2], u2[iu2], u3[iu2]]

        @NLconstraint(model, hs_defect_1(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_2(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_3(txu1..., txu2..., tf) == 0.0)
        # @NLconstraint(model, hs_defect_4(txu1..., txu2..., tf) == 0.0)
        # @NLconstraint(model, hs_defect_5(txu1..., txu2..., tf) == 0.0)
        # @NLconstraint(model, hs_defect_6(txu1..., txu2..., tf) == 0.0)
        @NLconstraint(model, hs_defect_7(txu1..., txu2..., tf) == 0.0)
    end
    return model, hs_defect, get_target_state
end