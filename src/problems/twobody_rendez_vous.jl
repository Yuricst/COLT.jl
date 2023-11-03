"""
Ready-to-go model for two-body spacecraft rendez-vous
"""


"""
Dynamics for two-body spacecraft rendez-vous
"""
function dynamics_twobody(t, states, u, p)
    x, y, z, _, _, _, mass = states  # unpack states
    r3 = (x^2 + y^2 + z^2)^(3/2)     # radius^3
    mu, c1, c2 = p                   # unpack parameters
    dx = [
        states[4],
        states[5],
        states[6],
        -mu/r3 * x + c1/mass*u[1],
        -mu/r3 * y + c1/mass*u[2],
        -mu/r3 * z + c1/mass*u[3],
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
    tf_bounds;
    mass_bounds=[0.5, 1.0],
    r_max::Real=3.0,
    v_max::Real=1.5,
    N::Int=40,
    solver=Ipopt.Optimizer,
    tf_free::Bool = false,
    max_mf::Bool = true,
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
        -r_max ≤ x[1:N_nodes]  ≤ r_max, (start = 1.0)
        -r_max ≤ y[1:N_nodes]  ≤ r_max, (start = 0.1)
        -r_max ≤ z[1:N_nodes]  ≤ r_max, (start = 0.0)
        -v_max ≤ vx[1:N_nodes] ≤ v_max, (start = 0.0)
        -v_max ≤ vy[1:N_nodes] ≤ v_max, (start = 1.1)
        -v_max ≤ vz[1:N_nodes] ≤ v_max, (start = 0.0)
        mass_bounds[1] ≤ mass[1:N_nodes] ≤ mass_bounds[2], (start = 1.0)
        # controls
        -1.0 ≤ u1[1:N_controls] ≤ 1.0, (start = 0.5)
        -1.0 ≤ u2[1:N_controls] ≤ 1.0, (start = 0.5)
        -1.0 ≤ u3[1:N_controls] ≤ 1.0, (start = 0.5)
    end);
    state_variables = [x,y,z,vx,vy,vz,mass]
    control_variables = [u1,u2,u3]

    # fix final time
    if tf_free == false
        fix(tf, tf_bounds[2]; force = true)
    end

    # fix node locations
    for i = 1:N_controls
        fix(t[i], (i-1)/(N_controls-1); force = true)   # fix the node locations 
    end

    # boundary constraints (initial)
    fix(x[1],    state0[1]; force = true)
    fix(y[1],    state0[2]; force = true)
    fix(z[1],    state0[3]; force = true)
    fix(vx[1],   state0[4]; force = true)
    fix(vy[1],   state0[5]; force = true)
    fix(vz[1],   state0[6]; force = true)
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
        get_target_state(_t) = keplerder_nostm_colt(mu, target_state, 0.0, _t, 1.e-12, 10)
        target_state = @expression(model, get_target_state(tf))
        @constraint(model, x[end]  - target_state[1] == 0.0)
        @constraint(model, y[end]  - target_state[2] == 0.0)
        @constraint(model, z[end]  - target_state[3] == 0.0)
        @constraint(model, vx[end] - target_state[4] == 0.0)
        @constraint(model, vy[end] - target_state[5] == 0.0)
        @constraint(model, vz[end] - target_state[6] == 0.0)
    end

    # path constraints
    for i = 1:N_controls
        @constraint(model, u1[i]^2 + u2[i]^2 + u3[i]^2 ≤ 1.0)
    end

    # objective is linear, to maximize final mass or minimize time
    if max_mf == true
        @objective(model, Max, mass[end]);
    else
        @objective(model, Min, tf);
    end
    return model, dynamics_twobody, state_variables, control_variables, [mu, c1, c2]
end