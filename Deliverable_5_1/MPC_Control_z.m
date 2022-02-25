classdef MPC_Control_z < MPC_Control
    properties
        A_bar, B_bar, C_bar % Augmented system for disturbance rejection
        L                   % Estimator gain for disturbance rejection
    end
    
    methods
        function mpc = MPC_Control_z(sys, Ts, H)
            mpc = mpc@MPC_Control(sys, Ts, H);
            
            [mpc.A_bar, mpc.B_bar, mpc.C_bar, mpc.L] = mpc.setup_estimator();
        end
        
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   d_est        - disturbance estimate
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps
            
            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.3)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar(1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % SET THE PROBLEM CONSTRAINTS con AND THE OBJECTIVE obj HERE
            
            
            A = mpc.A; B = mpc.B;
            Q = diag([1, 10]);
            R = eye(nu) * .1;
            % define the system 
            sys = LTISystem('A', A, 'B', B);
            % Input Pavg
            % z vz
            % Constraints on states
            umax = 80-56.667;
            umin = 50-56.667;
            xmax = [inf; inf];
          
            % cost function
            sys.x.penalty = QuadFunction(Q);
            sys.u.penalty = QuadFunction(R);
            
            % constraints
            sys.x.min = -xmax;
            sys.x.max = xmax;
            sys.u.max = umax;
            sys.u.min = umin;
            
            F = [eye(nx); -eye(nx)];
            f = [sys.x.max; -sys.x.min];
            M = [eye(nu); -eye(nu)];
            m = [sys.u.max; -sys.u.min];
            
            % extract LQR gain and terminal set
            Qf = sys.LQRPenalty.weight;
            % terminal set and cost
            sys.x.with('terminalPenalty');
            sys.x.terminalPenalty = QuadFunction(Qf);
            Bd = B;
            % Optimization constraints and objectives
            con = (X(:, 2) == A * X(:, 1) + Bd * U(:, 1) + B * d_est) + (M * U(:, 1) <= m);
            obj = (U(:, 1)-u_ref)' * R * (U(:, 1)-u_ref);
            
            
            for i = 2:N-1
                con = con + ((X(:,i+1)-x_ref)==A*(X(:,i)-x_ref)+B*(U(:,i)-u_ref));
                con = con + (F*X(:,i)<=f) + (M*U(:,i)<=m);
                obj = obj + (X(:,i)-x_ref)'*Q*(X(:,i)-x_ref) + (U(:,i)-u_ref)'*R*(U(:,i)-u_ref);
            end  
            obj = obj + (X(:,N)-x_ref)'*Qf*(X(:,N)-x_ref);
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref, d_est}, U(:,1));
        end
        
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            %   d_est  - disturbance estimate
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);
            
            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.3)
            ref = sdpvar;
            
            % Disturbance estimate (Ignore this before Part 5)
            d_est = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            A = mpc.A; B = mpc.B; C = mpc.C; D = mpc.D;
            Bd = B;
            Cd = 0;
            umax = 80-56.667;
            umin = 50-56.667;
            xmax = [inf; inf];
            F = [eye(nx); -eye(nx)];
            f = [xmax; xmax];
            M = [1; -1];
            m = [umax; -umin];
            con = [xs == A*xs + B*us + Bd*d_est, ref == C*xs + D*us + Cd*d_est, F*xs <= f, M*us <= m];
            
            obj = us'*us;   
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), {ref, d_est}, {xs, us});
        end
        
        
        % Compute augmented system and estimator gain for input disturbance rejection
        function [A_bar, B_bar, C_bar, L] = setup_estimator(mpc)
            
            %%% Design the matrices A_bar, B_bar, L, and C_bar
            %%% so that the estimate x_bar_next [ x_hat; disturbance_hat ]
            %%% converges to the correct state and constant input disturbance
            %%%   x_bar_next = A_bar * x_bar + B_bar * u + L * (C_bar * x_bar - y);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D
            [nx, nu] = size(mpc.B);
            ny = size(mpc.C,1);
            Bd = mpc.B;
            Cd = 0;
            A_bar = [mpc.A, Bd; zeros(nu,nx), eye(nu)]; %nx nx
            B_bar = [mpc.B; zeros(nu,nu)]; 
            C_bar = [mpc.C, Cd]; % ny nx
            poles = [0.4, 0.5, 0.6];
            L = -place(A_bar',C_bar', poles')';
            
            
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        
    end
end
