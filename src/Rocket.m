classdef Rocket
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Rocket class
    %
    % Structure (Click VIEW > Collapse All for better overview)
    %
    % methods
    % - System
    % - Simulation
    % - Visualization
    %
    % methods (hidden)
    % - System
    % - Simulation
    % - Visualization
    %
    % methods (static)
    % - Helpers
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    properties
        sys        % State, input, output names and units
        dimx = 12; % Length of state vector
        dimu = 4;  % Length of input vector
        Ts         % Sampling time
        
        anim_rate = 0; % Animation rate for visualization, default = 0x.
        
        mass   = 1.7;
    end
    properties (Constant)
        % Input limits
        ubu = [deg2rad([ 15;  15]); 80;  20]; % upper bound for input: [rad, rad, %, %]
        lbu = [deg2rad([-15; -15]); 50; -20]; % lower bound for input: [rad, rad, %, %]
        
        % State limits for (linear) control
        ubx = [deg2rad([Inf Inf Inf, 5 5 Inf]), Inf Inf Inf, Inf Inf Inf]';
        lbx = -[deg2rad([Inf Inf Inf, 5 5 Inf]), Inf Inf Inf, Inf Inf Inf]';
        
        % Indices in the state vector
        indx = struct('omega', 1:3, 'phi', 4:6, 'vel', 7:9, 'pos', 10:12);
        indu = struct('d1', 1, 'd2', 2, 'Pavg', 3, 'Pdiff', 4);
    end
    properties (Hidden, Constant)
        % Physical properties and constraints
        g      = 9.81;                                % gravitational acceleration
        J      = diag([0.0644; 0.0644; 0.0128]);      % inertia tensor
        inv_J  = inv(diag([0.0644; 0.0644; 0.0128]));
        
        % Propeller model constants
        thrust_coeff = [0; 0.03; 0];    % experimentally identified
        torque_coeff = -0.1040;         % experimentally identified
        r_F          = [0; 0; -0.215];  % thruster position in body frame
        
        % Plotting colors
        color = struct('meas', [0, 0.4470, 0.7410], 'ref', 'k');
    end
    
    methods
        %
        % Constructor
        %
        function obj = Rocket(Ts)
            
            obj.Ts = Ts;
            
            try % Test YALMIP installation
                sdpvar(2,1);
            catch
                error('Could not load YALMIP - check that it is installed properly and on the path.')
            end
            try % Test casadi installation
                import casadi.*
                casadi.SX.sym('x');
            catch
                error('Could not load casadi - check that it is installed properly and on the path.')
            end
            
            % Define system state, input, output names and inputs
            sys.StateName = { ...
                'wx', 'wy', 'wz', ...
                'alpha', 'beta', 'gamma',...
                'vx', 'vy', 'vz', ...
                'x','y','z'};
            sys.StateUnit = { ...
                'rad/s', 'rad/s', 'rad/s', ...
                'rad', 'rad', 'rad', ...
                'm/s', 'm/s', 'm/s', ...
                'm', 'm', 'm'};
            
            sys.InputName = {'d1', 'd2', 'Pavg', 'Pdiff'};
            sys.InputUnit = {'rad', 'rad', '%', '%'};
            
            sys.OutputName = sys.StateName;
            sys.OutputUnit = sys.StateUnit;
            
            obj.sys = sys;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % System
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Continuous rocket dynamics x_dot = f(x,u)
        % ----------------------------------------------------
        % x = [w; phi; v; p]
        % w     [rad/s] Angular velocity in body frame (b, {forward, left, up})
        % phi   [rad]   Euler angles expressing the attitude of body frame w.r.t. world frame (w, {East, North, Up})
        % v     [m/s]   Linear velocity in world frame
        % p     [m]     Position in world frame
        % ----------------------------------------------------
        % u = [d1; d2; P_avg; P_delta]
        % d1    [rad] Thrust vector deflection producing negativ angular velocity about body x axis
        % d2    [rad] Thrust vector deflection producing negativ angular velocity about body y axis
        % Pavg  [%]   Avg propeller throttle
        % Pdiff [%]   Difference propeller throttle producing negative angular velocity about body z axis
        %
        function [x_dot, y] = f(obj, x, u)
            
            if nargin < 3
                u = zeros(4, 1);
            end
            
            % Decompose state
            [w, phi, v, ~] = obj.parse_state(x);
            
            % Get transformation from body to world frame
            Twb = obj.eul2mat(phi);

            % Get force and moment in body frame from (differential) thrust
            [b_F, b_M] = obj.getForceAndMomentFromThrust(u);
            
            % State derivative components
            % Angular velocity in body frame
            w_dot = obj.inv_J * (b_M - cross(w, obj.J * w));
            
            % Attitude angles. Angular velocities in body frame -> rate of
            % change of Euler angles
            bet = phi(2); gam = phi(3);
            
            E_inv  = 1/cos(bet) * ...
                [cos(gam),         -sin(gam),          0;
                sin(gam)*cos(bet), cos(gam)*cos(bet),  0;
                -cos(gam)*sin(bet), sin(gam)*sin(bet), cos(bet)];
            
            phi_dot = E_inv * w;
            
            % Linear velocity and position in world frame   
            v_dot = Twb * b_F/obj.mass - [0; 0; obj.g];
            p_dot = v;
            
            % Compose state derivative
            x_dot = [w_dot; phi_dot; v_dot; p_dot];
            
            % Output thrust force and moment for debugging
            y = [b_F; b_M];
        end
        
        %
        % Get force and moment vector from thrust in body reference frame
        %
        function [b_F, b_M] = getForceAndMomentFromThrust(obj, u)
            
            % Thrust and torque magnitude from avg and diff throttle command
            thrust = obj.g * obj.thrust_coeff' * [u(3)^2; u(3); 1];
            torque = obj.torque_coeff * obj.J(3,3) * u(4);
            
            % Direction of thrust and differential torque vector in body frame (unit vector)
            b_eF = [sin(u(2)); -sin(u(1))*cos(u(2)); cos(u(1))*cos(u(2))];
            
            % Resulting force and torque vector in body frame
            b_F = thrust * b_eF;
            b_M = torque * b_eF + cross(obj.r_F, b_F);
        end
        
        %
        % Compute trim point (steady flight, x_dot = 0)
        %
        function [xs, us] = trim(obj)
            
            % Since we know here that a trim point f(xs, us) = 0 exists,
            % we can simply solve an box-constrained optimization problem
            % with Matlab's fmincon solver:
            %                   min f(xs,us).^2
            %
            % Given that we have nonlinear dynamics, we will find the local
            % minimum next to the initial guess.
            %
            % Decision variable y = [x; u]
            
            MSE = @(v) v'*v;
            objective_fcn = @(y) MSE( obj.f(y(1:12), y(13:16)) );
            
            % Set bounds and initial guess
            ubx_ = [Inf Inf Inf, deg2rad([180 89 180]), Inf Inf Inf, Inf Inf Inf]';
            lbx_ = -ubx_;
            
            uby = [ubx_; obj.ubu(:)];
            lby = [lbx_; obj.lbu(:)];
            
            y_guess = zeros(16,1);
            
            % Setup solver
            opt = optimoptions('fmincon', 'Algorithm','sqp');
            opt.Display = 'off';
            
            % Solve
            [y, fval, exitflag] = fmincon(objective_fcn, y_guess, ...
                [], [], [], [], lby, uby, [], opt);
            
            if exitflag < 0 || fval > 1e-3
                error('Could not find trim condition');
            end
            xs = y(1:12);
            us = y(13:16);
            
            % Clean up numerical inaccuracies
            xs(abs(xs) < 5e-2) = 0;
            us(abs(us) < 1e-3) = 0;
        end
        
        %
        % Return the linearization of the system around the
        % equilibrium point (xs, us)
        %
        function linSys = linearize(obj, xs, us)
            if nargin < 3
                fprintf('No equilibrium given... trimming\n');
                [xs, us] = obj.trim();
            end
            
            % Create casadi symbolic variables
            x = casadi.SX.sym('x', 12);
            u = casadi.SX.sym('u', 4);
            f = obj.f(x, u);
            
            % Create symbolic casadi function for automatic differentiation
            % of A = df/dx, B = df/du. Evaluate and densify.
            A = casadi.Function('A', {x,u}, {jacobian(f, x)});
            A = full(A(xs, us));
            B = casadi.Function('B', {x,u}, {jacobian(f, u)});
            B = full(B(xs, us));
            
            % Clean up numerical inaccuracies
            A(abs(A) < 1e-5) = 0;
            B(abs(B) < 1e-5) = 0;
            
            % Create state space representation
            linSys = ss(A, B, eye(12), zeros(12,4));
            linSys.UserData.xs = xs;
            linSys.UserData.us = us;
            
            linSys.StateName = obj.sys.StateName;
            linSys.StateUnit = obj.sys.StateUnit;
            
            linSys.InputName = obj.sys.InputName;
            linSys.InputUnit = obj.sys.InputUnit;
            
            linSys.OutputName = obj.sys.OutputName;
            linSys.OutputUnit = obj.sys.OutputUnit;
        end
        
        %
        % Decompose the system into four systems around a hovering
        % equilibrium
        %
        function [sys_x, sys_y, sys_z, sys_roll] = decompose(obj, linSys, xs, us)
            
            % Split into four seperate systems
            I = obj.indx;
            
            % [wy beta vx x], [d2]
            idx = [I.omega(2) I.phi(2) I.vel(1) I.pos(1)];
            idu = 2;
            idy = I.pos(1);
            
            sys_x = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys x');
            
            % [wx alpha vy y], [d1]
            idx = [I.omega(1) I.phi(1) I.vel(2) I.pos(2)];
            idu = 1;
            idy = I.pos(2);
            
            sys_y = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys y');
            
            % [vz z], [Pavg]
            idx = [I.vel(3) I.pos(3)];
            idu = 3;
            idy = I.pos(3);
            
            sys_z = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys z');
            
            % [wz gamma], [Pdiff]
            idx = [I.omega(3) I.phi(3)];
            idu = 4;
            idy = I.phi(3);
            
            sys_roll = obj.parse_system(linSys, xs, us, idx, idu, idy, 'sys roll');
        end
        
        %
        % Decompose the system into four systems around a hovering
        % equilibrium
        %
        function controller = merge_lin_controllers(obj, xs, us, mpc_x, mpc_y, mpc_z, mpc_roll)
            
            % Get state indices
            linSys = obj.linearize(xs, us);
            [sys_x, sys_y, sys_z, sys_roll] = obj.decompose(linSys, xs, us);
            
            Iu = zeros(4, 1);
            
            idx_x = sys_x.UserData.idx;
            idu_x = sys_x.UserData.idu;
            Iu_x  = Iu; Iu_x(idu_x) = 1;
            
            idx_y = sys_y.UserData.idx;
            idu_y = sys_y.UserData.idu;
            Iu_y  = Iu; Iu_y(idu_y) = 1;
            
            idx_z = sys_z.UserData.idx;
            idu_z = sys_z.UserData.idu;
            Iu_z  = Iu; Iu_z(idu_z) = 1;
            
            idx_r = sys_roll.UserData.idx;
            idu_r = sys_roll.UserData.idu;
            Iu_r  = Iu; Iu_r(idu_r) = 1;
            
            % If z_ is the state vector, 13:end = [], and mpc_z will be
            % evaluated with 0 disturbance
            controller = @(z_, ref_) us + ...
                Iu_x    .* mpc_x.get_u(  z_(idx_x) - xs(idx_x),              ref_(1) ) + ...
                Iu_y    .* mpc_y.get_u(  z_(idx_y) - xs(idx_y),              ref_(2) ) + ...
                Iu_z    .* mpc_z.get_u( [z_(idx_z) - xs(idx_z); z_(13:end)], ref_(3) ) + ...
                Iu_r .* mpc_roll.get_u(  z_(idx_r) - xs(idx_r),              ref_(4) );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Simulate any system in steps of Ts using an open-loop input
        % trajectory U, a single input u, or a control law @(x_, ref_)
        %
        function [T, X, U, Ref, Z_hat] = simulate(obj, f, x0, Tf, U, ref, esti)
            
            argno = struct('f',2, 'x0',3, 'Tf',4, 'U',5, 'ref',6, 'est',7);
            
            % Handle system dynamics
            if nargin < argno.f || (nargin >= argno.f && isempty(f))
                % Choose full nonlinear rocket dynamics
                f = @obj.f;
            end
            % Determine Ts, nx, simulation function
            if isa(f, 'function_handle')
                % Nonlinear model
                try x_dot = f(x0); catch, error('Length mismatch x0 / f(x,u)'); end
                nx = length(x_dot); clear dx;
                Ts_ = obj.Ts;
                sim_step = @(x_, u_) obj.simulate_step(f, x_, u_, obj.Ts);
                model_text = 'continuous-time nonlinear';
                
                if nargin < argno.U
                    U = zeros(obj.dimu, 1);
                    try f(x0, U); catch, error(['Cannot simulate autonomous system. ' ...
                            'Please provide a zeros input U of appropriate size.']); end
                end
                ubu_ = obj.ubu;
                lbu_ = obj.lbu;
                
            elseif isa(f, 'ss')
                % Linear model: Get size, setup sim_step
                lin_sys = f;
                nx = size(lin_sys.A, 1);
                if ~isfield(lin_sys.UserData, 'xs'), error('Specify xs in {state space model}.UserData.xs.'); end
                if ~isfield(lin_sys.UserData, 'us'), error('Specify us in {state space model}.UserData.us.'); end
                if ~isfield(lin_sys.UserData, 'idu'), error('Specify idu in {state space model}.UserData.idu.'); end
                
                ubu_ = obj.ubu(lin_sys.UserData.idu);
                lbu_ = obj.lbu(lin_sys.UserData.idu);
                
                if f.Ts == 0
                    % Continuous sys, discretize
                    Ts_ = obj.Ts;
                    fprintf('Discretize system ...')
                    sys_d = c2d(f, Ts_);
                else
                    sys_d = f;
                    Ts_ = sys_d.Ts;
                end
                sim_step = @(x_, u_) sys_d.A * x_ + sys_d.B * u_;
                model_text = 'discrete-time linear';
                %                 else
                %                     sys_c = f;
                %                     Ts_ = obj.Ts;
                %                     sim_step = @(x_, u_) obj.simulate_step( ...
                %                         @(x__, u__) sys_c.A * x__ + sys_c.B * u__, ...
                %                         x_, u_, obj.Ts);
                %                     model_text = 'continuous-time linear';
                %                 end
                if nargin < argno.U
                    % Autonomous system
                    U = zeros(lin_sys.B, 2);
                end
                
            else
                error(['First argument must be nonlinear dynamics ' ...
                    'function or state-space model.']);
            end
            
            % Determine simulation dimensions, allocate arrays
            T = 0:Ts_:Tf;
            nSteps = length(T);
            X = zeros(nx, nSteps + 1);
            
            % Set initial state
            X(:,1) = x0;
            
            % Handle input / control law
            if isnumeric(U)
                % U is static, run open-loop simulation -------------------
                loop_text = 'open';
                
                if size(U, 2) == 1
                    % U = u is a static input, replicate for simulation
                    U = repmat(U, 1, nSteps);
                else
                    % U is a static trajectory
                    if size(U, 2) ~= length(T)
                        error(['Static input trajectory must contain one ' ...
                            'column for each simulation step (' num2str(nSteps) ').']);
                    end
                end
                
                % Simulate open-loop along U trajectory
                fprintf(['  > Simulating ' model_text ' model in ' ...
                    loop_text '-loop to Tf = ' num2str(Tf) ':\n    ( ']); tic;
                for iStep = 1:nSteps
                    
                    % Print simulation time
                    if mod(T(iStep), 2) == 0, fprintf([num2str(T(iStep)) ' ']); end
                    
                    % Simulate next state from true state and input
                    X(:, iStep+1) = sim_step( X(:,iStep), U(:,iStep) );
                end
                if mod(T(iStep), 2) ~= 0, fprintf([num2str(T(iStep)) ' ']); end
                fprintf(')  Done (%.2fs)\n', toc);
                
                Ref   = [];
                Z_hat = [];
                
            elseif isa(U, 'function_handle')
                % U is a control law, run closed-loop simulation ----------
                loop_text = 'closed';
                ctrl_fcn = U;
                
                % Handle reference / reference function
                if nargin < argno.ref, error(['Using a control law requires ' ...
                        'to give a reference of appropriate size.']); end
                if isnumeric(ref) && size(ref, 2) == 1
                    % ref is a static reference
                    ref_fcn = @(t_, x_) ref;
                    
                elseif isa(ref, 'function_handle')
                    % ref is a function of time and state
                    ref_fcn = ref;
                    
                else
                    error(['Input a static reference or a reference ' ...
                        'function handle ref_fcn = @(t_, x).']);
                end
                
                % Handle estimator / estimation function
                if nargin >= argno.est && ~isempty(esti)
                    estimator_text = 'with estimator ';
                    est_fcn = esti.est_fcn;
                    nz = length(esti.z_hat0);
                    Z_hat = zeros(nz, nSteps);
                    Z_hat(:,1) = esti.z_hat0;
                else
                    % No estimator given, use true state
                    estimator_text = '';
                    est_fcn = @(z_hat_, u_, x_) x_;
                    nz = nx;
                    Z_hat = zeros(nz, nSteps);
                    Z_hat(:,1) = X(:,1);
                end
                
                % Determine ref dimension by evaluating once, allocate Ref
                ref_ = ref_fcn( 0, x0 );
                Ref = zeros(length(ref_), nSteps); clear ref_
                Ref(:,1) = ref_fcn( T(1), Z_hat(1:nx,1) );
                
                % Determine input dimension by evaluating once, allocate U
                u_ = ctrl_fcn( Z_hat(:,1), Ref(:,1) );
                U = zeros(length(u_), nSteps); clear u_
                
                % Simulate closed-loop using control law
                fprintf(['  > Simulating ' model_text ' model ' estimator_text ...
                    'in ' loop_text '-loop to Tf = ' num2str(Tf) ':\n    ( ']); tic;
                
                sim_success = true;
                for iStep = 1:nSteps
                    % Print simulation time
                    if mod(T(iStep), 2) == 0, fprintf([num2str(T(iStep)) ' ']); end
                    
                    % Compute reference ref from time t and state estimate x_hat
                    Ref(:,iStep) = ref_fcn( T(iStep), Z_hat(1:nx,iStep) );
                    % Compute input from augmented state estimate z_hat and reference ref
                    u = ctrl_fcn( Z_hat(:,iStep), Ref(:,iStep) );
                    U(:, iStep) = u;
                    if any(isnan(u))
                        sim_success = false;
                        break;
                    end
                        
                    % Simulate next true state from true state and input
                    X(:, iStep+1) = sim_step( X(:,iStep), U(:,iStep) );
                    
                    % Estimate augmented state from previous value, input u, and system measurement x
                    Z_hat(:,iStep+1) = est_fcn( Z_hat(:,iStep), U(:,iStep), X(:,iStep+1) );
                    
                end
                if size(Z_hat,1) == 13
                    figure
                    plot(Z_hat(end,:));
                    title('Estimate of disturbance');
                end
                
            
                if mod(T(iStep), 2) ~= 0, fprintf([num2str(T(iStep)) ' ']); end
                
                if sim_success
                    fprintf(')  Done (%.2fs)\n', toc);
                else
                    T = T(:,1:iStep);
                    X = X(:,1:iStep+1);
                    U = U(:,1:iStep);
                    Ref = Ref(:,1:iStep);
                    Z_hat = Z_hat(:,1:iStep+1);
                    fprintf(')  Abort (%.2fs)\n', toc);
                    fprintf('Problem when computing control from x = [');
                    fprintf(' %g', X(:,iStep));
                    fprintf(' ]\n');
                end
                
                % Remove last state for consistent size with U
                Z_hat(:,end) = [];
            else
                error(['Input a static input or a control ' ...
                    'function handle ctrl_fcn = @(x_, ref_).']);
            end
            
            % Remove last state for consistent size with U
            X(:,end) = [];
            
            % If linear system, return absolute state and input
            if isa(f, 'ss')
                X = X + lin_sys.UserData.xs;
                U = U + lin_sys.UserData.us;
            end
            
            % (Assume we know about the system limits when ?bu_ length matches U)
            % Warn if input constraints have been driven above physical
            % system limits (plus numerical tolerance)
            if length(lbu_) == size(U, 1)
                num_tol = [deg2rad([0.1, 0.1]) 0.1 0.1]';
                lb_exceeded = any(U < (lbu_ - num_tol), 2); ub_exceeded = any((ubu_ + num_tol) < U, 2);
                warning('off','backtrace')
                if any([ub_exceeded; lb_exceeded])
                    for i = 1:length(ub_exceeded)
                        if ub_exceeded(i)
                            warning(['Upper physical limit of input ' num2str(i) ' has been exceeded during the simulation.']);
                        end
                        if lb_exceeded(i)
                            warning(['Lower physical limit of input ' num2str(i) ' has been exceeded during the simulation.']);
                        end
                    end
                end
                warning('on','backtrace')
            end
        end
        
        %
        % Simulate the nonlinear model in steps of obj.Ts using an open-loop
        % input trajectory U, a single input u, or a control law @(x_, ref_)
        %
        function [T, X, U, Ref, Z_hat] = simulate_f(obj, x0, Tf, U, ref, esti)
            
            if nargin >=6
                [T, X, U, Ref, Z_hat] = obj.simulate([], x0, Tf, U, ref, esti);
            elseif nargin >= 5
                [T, X, U, Ref, Z_hat] = obj.simulate([], x0, Tf, U, ref);
            elseif nargin >= 4
                [T, X, U, Ref, Z_hat] = obj.simulate([], x0, Tf, U);
            else
                [T, X, U, Ref, Z_hat] = obj.simulate([], x0, Tf);
            end
            
        end
        
        %
        % Simulate the nonlinear model in steps of obj.Ts using an open-loop
        % input trajectory U, a single input u, or a control law @(x_, ref_).
        % Use hardcoded observer in z direction.
        %
        function [T, X, U, Ref, Z_hat] = simulate_f_est_z(obj, x0, Tf, U, ref, mpc_z, sys_z)
            
            esti.est_fcn = @(z_hat_, u_, x_) estimator_z(z_hat_, u_, x_, mpc_z, sys_z);
            esti.z_hat0 = [x0; 0]; % Initial guess for estimation
            
            [T, X, U, Ref, Z_hat] = obj.simulate([], x0, Tf, U, ref, esti);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Visualization
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Visualize the trajectory (without plots)
        %
        function ph = vis(obj, T, X, U, ax)
            
            if nargin < 5
                % If no plot handles are provided, create new plot
                ph = obj.create_pose_axes();
            else
                ph = obj.create_pose_axes(ax);
            end
            ax = ph.ax_pose;
            %cla(ax);
            
            if length(T) < 2
                sim_dt = 0.05;
            else
                sim_dt = T(2) - T(1);
            end
            nT = length(T);
            
            if size(U, 2) < 2
                % If U is a vector, replicate it to match state trajectory
                U = repmat(U, 1, nT);
            end
            
            if obj.anim_rate > 0
                % Animate trajectory
                fprintf('  > Visualizing (%.1fx) ... ', obj.anim_rate);
                
                fps_max = 20; % Max frames per second to display
                vis_dt_ref = max(1/fps_max, sim_dt/obj.anim_rate);
                iStep = ceil(vis_dt_ref / (sim_dt/obj.anim_rate));
                vis_dt = iStep * sim_dt/obj.anim_rate;
                
                % Create rocket object
                [ph_rocket, ph_thrust] = obj.create_rocket_transformobj(ax);
                
                % Transform it along trajectory / draw thrust
                for iT = 1:iStep:length(T)
                    tic;
                    % Plot missed position points
                    if iT > 1
                        obj.visualize_pos( X(obj.indx.pos, iT-iStep:iT), ax);
                    end
                    % Update properties of rocket object
                    obj.visualize_point( X(:,iT), U(:,iT), ax, ph_rocket, ph_thrust );
                    % Pause by remaining time after ploting point
                    time_lapsed = toc;
                    pausetime = vis_dt - time_lapsed;
                    pause(pausetime);
                end
                fprintf('Done\n')
            else
                % Plot full position trajectory
                [~, ~, ~, pos] = obj.parse_state(X);
                obj.visualize_pos(pos, ax);
                
                % Plot pose at discrete points in time (in any case initial and final state)
                snapshots_per_sim_sec = 1.1;    % (Min) number of snapshots per simulation second
                iStep = floor(1/(snapshots_per_sim_sec * sim_dt));
                iStart = mod(nT, iStep);
                if iStart == 0, iStart = iStep; end
                
                % Plot initial pose, then continue at iStart at iSteps to
                % end up exactly on the last time step.
                obj.visualize_point( X(:,1), U(:,1), ax );
                for iT = iStart:iStep:nT
                    obj.visualize_point( X(:,iT), U(:,iT), ax );
                end
            end
            ph.ax_pose = ax;
        end
        
        %
        % Plot & visualize full trajectory
        %
        function ph = plotvis(obj, T, X, U, Ref)
            
            X_ref = nan(size(X));
            if nargin >= 5
                if size(Ref, 2) == 1
                    Ref = repmat(Ref, 1, length(T));
                end
                
                if size(Ref, 1) == 12
                    % Full-state reference is already given
                    X_ref = Ref;
                elseif size(Ref, 1) == 4
                    % Got ref4
                    X_ref = nan(size(X));
                    id_ref = [obj.indx.pos, obj.indx.phi(3)];
                    X_ref(id_ref, :) = Ref;
                else
                    warning('Cannot map reference to states. Provide a reference with 4 or 12 columns.');
                end
            end
            
            % Create axes for full trajectory
            nxu = 12 + 4;
            nColVis = 2; % Number of columns for pose visualization
            
            targetSize = [1400, 730];
            position = obj.get_position_on_screen(targetSize);
            ph.fig = figure('Position', position);

            nrow = 5;
            ncol = 3 + nColVis;
            
            ax = gobjects(nxu, 1);
            
            % Map input/states to subplot axes
            u_id = [1 1 2 3];
            x_id = [
                1*ncol+[1 1 3] % wy wx wz
                2*ncol+[1 1 3] % beta alpha gamma
                3*ncol+[1 1 2] % vx vy vz
                4*ncol+[1 1 2] % x y z
                ]';
            x_id = x_id(:)';
            ux_id = [u_id x_id];
            
            % Create subplot axes
            for ixu = 1:nxu
                ax(ux_id(ixu)) = subplot(nrow, ncol, ux_id(ixu));
                ax_ = ax(ux_id(ixu));
                hold(ax_, 'on');
                grid(ax_, 'on');
                axis(ax_, 'tight');
            end, clear ixu ax_
            
            % Write title and time
            title(ax(1), '~ sys x (blue), y (red)');
            title(ax(2), '~ sys z');
            title(ax(3), '~ sys roll');
            
            for i = [2*ncol+3, 4*ncol + [1 2]]
                xlabel(ax(i), 'Time [s]');
            end, clear i
            
            % Plot data
            bx = [
                -[deg2rad([Inf Inf Inf, 5 Inf Inf]), Inf Inf Inf, Inf Inf Inf];
                [deg2rad([Inf Inf Inf, 5 Inf Inf]), Inf Inf Inf, Inf Inf Inf]
                ]';
            bu = [obj.lbu(:), obj.ubu(:)];
            
            obj.plot_into_axes(ax, ux_id, obj.sys, T, X, U, bx, bu, X_ref);
            
            % Visualize data
            ax_pose = subplot(nrow, ncol, [ncol - nColVis + 1, nrow*ncol]);
            obj.vis(T, X, U, ax_pose);
            
            % Plot reference if not constant
            pos_ref = X_ref(obj.indx.pos, :);
            if size(unique(pos_ref', 'rows'), 1) > 1
                plot3(ax_pose, pos_ref(1,:), pos_ref(2,:), pos_ref(3,:), ...
                    'Color', obj.color.ref, 'LineStyle', '--', ...
                    'LineWidth', 1);
            end
            
            ph.ax_ux = ax;
            ph.ax_pose = ax_pose;
        end
        
        %
        % Plot and visualize sub-system trajectory (fill up other states
        % with (xs, us) for 3D visualization
        %
        function ph = plotvis_sub(obj, T, subX, subU, subSys, xs, us, ref)
            
            if nargin < 8
                ref = NaN;
            end
            
            % Create axes for two column figure
            nxu = size(subX, 1) + size(subU, 1);
            nColVis = 1; % Number of columns for pose visualization
            
            targetSize = [1120, 420];
            position = obj.get_position_on_screen(targetSize);
            ph.fig = figure('Position', position);

            nrow = nxu;
            ncol = 1 + nColVis;
            
            ax = gobjects(nxu, 1);
            
            % Map inputs/states to subplot axes
            ux_id = 1:ncol:ncol*nxu;
            
            % Create subplot axes
            for ixu = 1:nxu
                ax(ux_id(ixu)) = subplot(nrow, ncol, ux_id(ixu));
                ax_ = ax(ux_id(ixu));
                hold(ax_, 'on');
                grid(ax_, 'on');
                axis(ax_, 'tight');
            end, clear ixu ax_
            
            % Write title and time
            title(ax(1), subSys.Name);
            xlabel(ax(max(ux_id)), 'Time [s]');
            
            % Create full data
            if size(subX, 1) < 13
                idx = subSys.UserData.idx;
                X = repmat(xs, 1, length(T));
                X(idx,:) = subX;
            end
            if size(subU, 1) < 4
                idu = subSys.UserData.idu;
                U = repmat(us, 1, length(T));
                U(idu,:) = subU;
            end
            
            % Plot data
            bx = [obj.lbx(:), obj.ubx(:)];
            sub_bx = bx(idx,:);
            
            bu = [obj.lbu(:), obj.ubu(:)];
            sub_bu = bu(idu,:);
            
            % Create ref trajectory with ref in ref row and Nan otherwise
            subX_ref = nan(size(subX));
            subX_ref(subSys.UserData.idx == subSys.UserData.idy,:) = repmat(ref, 1, length(T));
            
            obj.plot_into_axes(ax, ux_id, subSys, T, subX, subU, sub_bx, sub_bu, subX_ref);
            
            % Visualize data
            ax_pose = subplot(nrow, ncol, [ncol - nColVis + 1, nrow*ncol]);
            obj.vis(T, X, U, ax_pose);
            
            ph.ax_xu = ax;
            ph.ax_pose = ax_pose;
        end
        
    end
    
    methods (Hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % System
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Continuous rocket dynamics x_dot = f(x, u) along 1D z direction
        %
        function x_dot = f_z(obj, x, u)
            
            if nargin < 3
                u = 1;
            end
            
            % Decompose state
            vz = x(1);
            
            thrust = obj.g * obj.thrust_coeff' * [u^2; u; 1];
            
            % State derivative components
            % Vertical speed and height
            vz_dot = thrust/obj.mass - obj.g;
            z_dot = vz;
            
            % Compose state derivative
            x_dot = [vz_dot; z_dot];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Simulation (hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Simulate the system xdot = f(x,u) from x0 forward for dt seconds
        %
        function xp = simulate_step(~, f, x0, u, dt)
            % Integrate forward to next time step
            [~, xout] = ode45( @(t_, x_) f(x_, u), [0, dt], x0 );
            xp = xout(end,:)';
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Visualization (hidden)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a graphical representation of the rocket
        %
        function [h_rocket, h_thrust] = create_rocket_transformobj(obj, ax)
            
            % Geometric properties
            geom.h = 0.6;  % [m]
            geom.r = 0.15; % [m]
            
            % Create hgtransform object for rocket
            if nargin < 2
                ph = obj.create_pose_axes();
                ax = ph.ax_pose;
                h_rocket = hgtransform('Parent', ax);
            else
                h_rocket = hgtransform('Parent', ax);
            end
            
            % Plot body axes ----------------------------------------------
            axes_scale = 0.3;
            
            axes_length = axes_scale * 1;
            axes_width  = axes_scale * 20;
            
            Xax = [zeros(3, 1), axes_length * [1;0;0]];
            Yax = [zeros(3, 1), axes_length * [0;1;0]];
            Zax = [zeros(3, 1), axes_length * [0;0;1]];
            line(ax, Xax(1,:), Yax(1,:), Zax(1,:), 'Color', 'r', 'LineWidth', axes_width, 'Parent', h_rocket);
            line(ax, Xax(2,:), Yax(2,:), Zax(2,:), 'Color', 'g', 'LineWidth', axes_width, 'Parent', h_rocket);
            line(ax, Xax(3,:), Yax(3,:), Zax(3,:), 'Color', 'b', 'LineWidth', axes_width, 'Parent', h_rocket);
            
            % Plot center of gravity (foreground) -------------------------
            plot3(ax, 0, 0, 0, ...
                '.', 'Color', 'k', 'MarkerSize', 50, 'Parent', h_rocket);
            
            % Plot zero thrust deflection direction -----------------------
            %X = obj.r_F(1) + [0 0];
            %Y = obj.r_F(2) + [0 0];
            %Z = obj.r_F(3) + [0 -1];
            %line(ax, X, Y, Z, 'LineStyle', '--', 'Color', [0 0 0 0.1], 'LineWidth', 4, 'Parent', h_rocket);
            
            % Plot thrust vector ------------------------------------------
            % ... with 'particle' markers
            h_thrust(1) = plot3(ax, obj.r_F(1), obj.r_F(2), obj.r_F(3), ...
                '*', 'Color', '#0072BD', 'MarkerSize', 8, 'LineWidth', 1, 'Parent', h_rocket);
            
            % ... with line
            h_thrust(2) = line(ax, obj.r_F(1), obj.r_F(2), obj.r_F(3), 'Color', '#0072BD', 'LineWidth', 4, 'Parent', h_rocket);
            
            % Outlet position
            plot3(ax, obj.r_F(1), obj.r_F(2), obj.r_F(3), ...
                '*', 'Color', '#D95319', 'MarkerSize', 12, 'LineWidth', 2, 'Parent', h_rocket);
            
            % Plot cylindric body shape -----------------------------------
            [X,Y,Z] = cylinder(geom.r, 13);
            Z = geom.h * Z + obj.r_F(3);
            surf(ax, X, Y, Z, ...
                'EdgeColor','k', 'EdgeAlpha', 0.2, ...
                'FaceColor','k', 'FaceAlpha', 0.1, ...
                'Parent', h_rocket);
        end
        
        %
        % Plot center of gravity point(s) into ax
        %
        function visualize_pos(obj, pos, ax)
            plot3(ax, pos(1,:), pos(2,:), pos(3,:), '.', ...
                'Color', obj.color.meas, ...
                'MarkerSize', 10);
        end
        
        %
        % Visualize the rocket at a given state and input
        %
        function ph = visualize_point(obj, x, u, ax, h_rocket, h_thrust)
            
            if nargin < 4
                % If no plot handles are provided, create new plot
                ph.fig = figure;
                ax = axes('Parent', ph.fig);
                ph.ax_pose = ax;
                axis(ax, 'equal');
                view(ax, 3)
                hold(ax, 'on')
                grid(ax, 'on')
            else
                ph.fig = ax.Parent;
                ph.ax = ax;
            end
            
            if nargin < 6
                % If no rocket graph obj available, create it once
                [h_rocket, h_thrust] = obj.create_rocket_transformobj(ax);
            end
            
            [~, phi, ~, pos] = obj.parse_state(x);
            
            % Rotation from body to interial frame
            Twb = obj.eul2mat(phi);
            h_rocket.Matrix(1:3, :) = [Twb, pos];
            
            % Plot thrust vector ------------------------------------------
            thrust_max = 29.43; % = obj.g * obj.thrust_coeff' * [100^2; 100; 1];
            thrust_scale = 1; % ENTER visualized vector length at max thrust
            
            % Thrust vector, scale for visualization
            [b_F, ~] = obj.getForceAndMomentFromThrust(u);
            Fvis = thrust_scale/thrust_max * -b_F;
            
            % ... with markers
            % %             thrust_resolution = 15;
            % %             Fx = linspace(0, Fvis(1), thrust_resolution);
            % %             Fy = linspace(0, Fvis(2), thrust_resolution);
            % %             Fz = linspace(0, Fvis(3), thrust_resolution);
            % %             h_thrust = h_thrust(1);
            % ... with line (more efficient)
            Fx = [0, Fvis(1)];
            Fy = [0, Fvis(2)];
            Fz = [0, Fvis(3)];
            h_thrust = h_thrust(2);
            
            h_thrust.XData = obj.r_F(1) + Fx;
            h_thrust.YData = obj.r_F(2) + Fy;
            h_thrust.ZData = obj.r_F(3) + Fz;
            
            % Plot center of gravity trace --------------------------------
            plot3(ax, pos(1), pos(2), pos(3), '.', ...
                'Color', obj.color.meas, ...
                'MarkerSize', 10);
        end
        
        %
        % Plot content into predefined axes
        %
        function plot_into_axes(obj, ax, ux_id, sys, T, X, U, bx, bu, X_ref)
            
            lw = 1;
            
            nx = size(X, 1);
            nu = size(U, 1);
            
            % Change all RAD content to DEG
            function [unit, var, bounds] = RAD_to_DEG(unit, var, bounds)
                for i = 1:length(unit)
                    if contains(unit{i}, 'rad')
                        unit{i} = strrep(unit{i}, 'rad', 'deg');
                        var(i,:) = rad2deg(var(i,:));
                        bounds(i,:)   = rad2deg(bounds(i,:));
                    end
                end
            end
            
            [sys.InputUnit, U, bu] = RAD_to_DEG(sys.InputUnit, U, bu);
            [~, X_ref, ~] = RAD_to_DEG(sys.StateUnit, X_ref, bx);
            [sys.StateUnit, X, bx] = RAD_to_DEG(sys.StateUnit, X, bx);
            
            % Plot constraints and axis labels
            function plot_constraints_and_yLabels(ax_, lub, axisName, axisUnit)
                
                line(ax_, repmat([T(1) T(end)], 2, 1)', repmat(lub, 2, 1), 'Color', 'k', 'LineStyle', '-.', 'LineWidth', lw);
                
                YLabel = ax_.YLabel.String;
                YLabeli = {[axisName ' [' axisUnit ']']};
                
                if isempty(YLabel)
                    YLabel = cell(YLabeli);
                else
                    YLabel = [YLabel, YLabeli];
                end
                
                ylabel(ax_, YLabel);
            end
            
            % Plot input and state constraints
            for iu = 1:nu, plot_constraints_and_yLabels(ax(ux_id(iu)),    bu(iu,:), sys.InputName{iu}, sys.InputUnit{iu}); end, clear iu
            for ix = 1:nx, plot_constraints_and_yLabels(ax(ux_id(nu+ix)), bx(ix,:), sys.StateName{ix}, sys.StateUnit{ix}); end, clear ix
            
            % Plot input
            for iu = 1:nu, plot( ax(ux_id(iu)),  T, U(iu, :), 'LineWidth', lw ); end, clear iu
            % Plot state
            for ix = 1:nx, plot( ax(ux_id(nu+ix)), T, X(ix, :), 'LineWidth', lw ); end, clear ix
            % Plot reference
            for ix = 1:nx, plot( ax(ux_id(nu+ix)), T, X_ref(ix,:), 'Color', obj.color.ref, 'LineStyle', '--', 'LineWidth', lw ); end, clear ix
        end
        
        %
        % Get safe outer figure size for specific screen
        %
        function Position = get_position_on_screen(~, targetSize, max_screen_usage)

            if nargin < 3
                max_screen_usage = 0.9;
            end

            figRatio = targetSize(1) / targetSize(2);

            screenSize = get(0, 'ScreenSize'); screenSize = screenSize(3:4);
            availSize = max_screen_usage * screenSize;
            horClip = [availSize(1), availSize(1) / figRatio];
            verClip = [figRatio * availSize(2), availSize(2)];
            minClip = min(horClip, verClip);
            figSize = round( min(minClip, targetSize) );
            
            % Determine figure offset (window/menu bar)
            f = figure('Visible', 'off');
            sizeOffset = f.OuterPosition(3:4) - f.Position(3:4);
            close(f);

            outerPos = screenSize - figSize - sizeOffset - [0 50]; % 50 height safety for Windows
            Position = [outerPos, figSize];
        end
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Helpers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Split the state into its parts, in: State trajectory (nx x time)
        %
        function [omega, phi, vel, pos] = parse_state(X)
            if nargout >= 1, omega = X(Rocket.indx.omega, :); end
            if nargout >= 2, phi   = X(Rocket.indx.phi, :); end
            if nargout >= 3, vel   = X(Rocket.indx.vel,   :); end
            if nargout >= 4, pos   = X(Rocket.indx.pos,   :); end
        end
        
        %
        % Obtain transformation matrix from attitude angles
        %
        function Twb = eul2mat(eul)
            alp = eul(1); % alpha
            bet = eul(2); % beta
            gam = eul(3); % gamma
            
            % T1 is elementary rotation about x axis
            function T = T1(a)
                T = [1 0 0;
                    0 cos(a) sin(a);
                    0 -sin(a) cos(a)];
            end
            % T2 is elementary rotation about y axis
            function T = T2(a)
                T = [cos(a) 0 -sin(a);
                    0 1 0;
                    sin(a) 0 cos(a)];
            end
            % T3 is elementary rotation about z axis
            function T = T3(a)
                T = [cos(a) sin(a) 0;
                    -sin(a) cos(a) 0;
                    0 0 1];
            end
            % %             % Detailed derivation
            % %
            % %             % Going from world to body frame:
            % %             % 1. Rotate by alpha about world x axis:
            % %             T_bw1 = T1(alp);
            % %             % 2. Rotate by beta about resulting y axis;
            % %             T_bw2 = T2(bet) * T_bw1;
            % %             % 3. Rotate by gamma about resluting z axis;
            % %             T_bw = T3(gam) * T_bw2;
            % %
            % %             % In one step:
            % %             T_bw = T3(gam) * T2(bet) * T1(alp);
            % %
            % %             % T_bw transforms a w_vect in body frame (b_vect):
            % %             % b_vect = T_bw * w_vect
            % %             % We need the opposite, so we transpose the matrix:
            % %             T_wb = T_bw';
            % %
            % %             % Alternatively, we can directly construct this matrix by going
            % %             % backwards from body to world frame (same logic, inverse
            % %             % angles):
            Twb = T1(-alp) * T2(-bet) * T3(-gam);
        end
        
        
        
        %
        % Create pose axes (3D, equal, hold, grid, labels)
        %
        function ph = create_pose_axes(figNo_or_ax)
            
            if nargin < 1 || (nargin >= 1 && ~isgraphics(figNo_or_ax))
                % Figure number is not specified, create figure and axes
                ph.fig = figure;
                ax = axes('Parent', ph.fig);
            else
                % Use given axes
                ax = figNo_or_ax;
                ph.fig = ax.Parent;
            end
            
            axis(ax, 'equal');
            hold(ax, 'on');
            grid(ax, 'on');
            view(ax, 3);
            
            title(ax, 'Pose');
            xlabel(ax, 'x [m]')
            ylabel(ax, 'y [m]')
            zlabel(ax, 'z [m]')
            ph.ax_pose = ax;
        end
        
        %
        % Create sub-system from indices
        %
        function sub_sys = parse_system(sys, xs, us, idx, idu, idy, name)
            
            [A, B, C, ~] = ssdata(sys);
            
            sub_sys = ss(A(idx,idx), B(idx,idu), C(idy,idx), 0);
            
            sub_sys.UserData.idx = idx;
            sub_sys.UserData.idu = idu;
            sub_sys.UserData.idy = idy;
            
            sub_sys.UserData.xs = xs(idx);
            sub_sys.UserData.us = us(idu);
            
            sub_sys.Name = name;
            sub_sys.StateName = sys.StateName(idx);
            sub_sys.StateUnit = sys.StateUnit(idx);
            
            sub_sys.InputName = sys.InputName(idu);
            sub_sys.InputUnit = sys.InputUnit(idu);
            
            sub_sys.OutputName = sys.OutputName(idy);
            sub_sys.OutputUnit = sys.OutputUnit(idy);
        end
        
        %
        % Trace out an MPC in ref_time seconds
        %
        function Ref4 = MPC_ref(t, ref_time, roll_max, tilt)
            
            if nargin < 3
                roll_max = deg2rad(15);
            end
                
            if nargin < 4
                tilt = true;
            end
                
            % Coordinates (x, y, heading)
            coords = [ ...
                0 0  45; 0 2  45; 1 1  45; 2 2  45; 2 0 -90;          % 'M'
                3 0 -90; 3 2 -90; 4 2 -90; 4 1 -90; 3 1 -90; 3 0 -90; % 'P'
                7 0 0; 5 0 90; 5 2 45; 7 2 0];                        % 'C'
            coords(:,1:2) = coords(:,1:2) / 2;
            coords(:,3) = coords(:,3) * rad2deg(roll_max)/90;
            
            nCoords = size(coords, 1);
            
            % Break the path into legs, compute their end times such that
            % the end of the path is reached at ref_time
            legs = coords(2:end,1:2) - coords(1:end-1,1:2);
            distances = vecnorm(legs, 2, 2);
            leg_endtimes = [0; ref_time * cumsum(distances) / sum(distances)]';
            
            % Find target index for each time point
            target_id = sum(t(:) > leg_endtimes, 2) + 1;
            % Limit target_id to final point
            target_id = min(nCoords, target_id);
            
            % Return target coordinates for each time point
            XZ = coords(target_id,1:2);
            Roll = deg2rad(coords(target_id,3));
            Ref4 = [XZ(:,1), 0, XZ(:,2), Roll]; % 4D X 0 Z roll
            
            % Rotate about x axis
            if tilt
                alpha = deg2rad(20);
                XZ = ([cos(alpha) -sin(alpha); sin(alpha) cos(alpha)] * XZ')';
                
                % Rotate about z axis
                gamma = deg2rad(-30);
                Ref4 = [cos(gamma) * XZ(:,1), sin(gamma) * XZ(:,1), XZ(:,2), Roll]; % 4D X 0 Z roll
            end
        end
    end
    
end