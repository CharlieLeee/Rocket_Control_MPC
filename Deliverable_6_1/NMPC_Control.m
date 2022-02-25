function opti_eval = NMPC_Control(rocket, H)

    import casadi.*
    opti = casadi.Opti(); % Optimization problem

    N = ceil(H/rocket.Ts); % MPC horizon
    nx = 12; % Number of states
    nu = 4;  % Number of inputs

    % Decision variables (symbolic)
    X_sym = opti.variable(nx, N); % state trajectory
    U_sym = opti.variable(nu, N-1);   % control trajectory)

    % Parameters (symbolic)
    x0_sym  = opti.parameter(nx, 1);  % initial state
    ref_sym = opti.parameter(4, 1);   % target position

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
    
    f_discrete = @(x,u) f_discrete_(rocket, x, u);
    Q = diag([40,40,1,10,10,20,1,1,1,20,20,20]);
    %Q = 10*eye(nx);
    R = 0.1*eye(nu);
    
    % steady state
    [xs_t, us_t] = rocket.trim(); % ref pose is vertical to ground
    xs = [xs_t(1:5)',ref_sym(4)', xs_t(7:9)', ref_sym(1:3)']';
    us = [us_t(1:end)']';
    % Terminal Cost and Terminal Set
    lin_sys = rocket.linearize(xs_t, us_t);
    lin_sys_d = c2d(lin_sys, rocket.Ts);    
    sys = LTISystem('A', lin_sys_d.A, 'B', lin_sys_d.B);
    sys.x.penalty = QuadFunction(Q); 
    sys.u.penalty = QuadFunction(R);
    sys.u.min = rocket.lbu; sys.u.max = rocket.ubu;
    sys.x.min = rocket.lbx; sys.x.max = rocket.ubx;
    Qf = sys.LQRPenalty.weight;
    Xf = sys.LQRSet % Empty
    
    %[~,Qf] = dlqr(lin_sys_d.A, lin_sys_d.B, Q, R);
    %Qf = dlyap(lin_sys_d.A,Q);
    
    % Terminal Condition and Cost
    obj = (xs - X_sym(:,N))' * Qf * (xs - X_sym(:,N));
    %opti.subject_to(Xf.A * X_sym(:,N) <= Xf.b);
        
    % 
    for k = 1:N-1
        obj = obj + (xs - X_sym(:,k))' * Q * (xs - X_sym(:,k))...
            + (us - U_sym(:, k))' * R * (us - U_sym(:, k));
        opti.subject_to(X_sym(:, k+1) == f_discrete(X_sym(:, k), U_sym(:,k)));
    end
    
    opti.minimize(obj);
    
    % constraints
    %opti.subject_to(-deg2rad(85) <= X_sym(4) <= deg2rad(85));
    opti.subject_to(-deg2rad(85) <= X_sym(5) <= deg2rad(85));
    opti.subject_to(-deg2rad(15) <= U_sym(1) <= deg2rad(15));
    opti.subject_to(-deg2rad(15) <= U_sym(2) <= deg2rad(15));
    opti.subject_to(50 <= U_sym(3) <= 80);
    opti.subject_to(-20 <= U_sym(4) <= 20);
    
    % boundary
    opti.subject_to(X_sym(:,1) == x0_sym);
        

    % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ---- Setup solver ------
    ops = struct('ipopt', struct('print_level', 0, 'tol', 1e-3), 'print_time', false);
    opti.solver('ipopt', ops);

    % Create function to solve and evaluate opti
    opti_eval = @(x0_, ref_) solve(x0_, ref_, opti, x0_sym, ref_sym, U_sym);
end

function u = solve(x0, ref, opti, x0_sym, ref_sym, U_sym)

    % ---- Set the initial state and reference ----
    opti.set_value(x0_sym, x0);
    opti.set_value(ref_sym, ref);

    % ---- Solve the optimization problem ----
    sol = opti.solve();
    assert(sol.stats.success == 1, 'Error computing optimal input');

    u = opti.value(U_sym(:,1));

    % Use the current solution to speed up the next optimization
    opti.set_initial(sol.value_variables());
    opti.set_initial(opti.lam_g, sol.value(opti.lam_g));
end

function xn = f_discrete_(rocket, x, u)
    %[x_dot, y] = rocket.f(x,u);
    %xn = rocket.Ts * x_dot + x;
    f = @(x,u) rocket.f(x,u);
    xn = RK4(x, u, rocket.Ts, f);
end

