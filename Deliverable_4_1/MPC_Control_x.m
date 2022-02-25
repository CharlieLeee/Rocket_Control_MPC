classdef MPC_Control_x < MPC_Control
    
    methods
        % Design a YALMIP optimizer object that takes a steady-state state
        % and input (xs, us) and returns a control input
        function ctrl_opti = setup_controller(mpc, Ts, H)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   X(:,1)       - initial state (estimate)
            %   x_ref, u_ref - reference state/input
            % OUTPUTS
            %   U(:,1)       - input to apply to the system
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            N = ceil(H/Ts); % Horizon steps

            [nx, nu] = size(mpc.B);
            
            % Targets (Ignore this before Todo 3.2)
            x_ref = sdpvar(nx, 1);
            u_ref = sdpvar(nu, 1);
            
            % Predicted state and input trajectories
            X = sdpvar(nx, N);
            U = sdpvar(nu, N-1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            
            % NOTE: The matrices mpc.A, mpc.B, mpc.C and mpc.D are
            %       the DISCRETE-TIME MODEL of your system
            
            % SET THE PROBLEM CONSTRAINTS con AND THE obj obj HERE

            x_indexes = [2,5,7,10];%wy, beta, vx, x
            u_index = 2;
            Q = diag([40,10,1,20]);
            R = 0.1;

            global rocket
            % Constraints
            % u in U = { u | Mu <= m }
            M = [1;-1]; m = [rocket.ubu(u_index); -rocket.lbu(u_index)];
            % x in X = { x | Fx <= f }
            F = zeros(nx*2, nx);
            f = zeros(nx*2, 1);
            for i = 1:nx
                F(2*i - 1, i) = 1;
                F(2*i , i) = -1;
                f(2*i-1, 1) = rocket.ubx(x_indexes(i));
                f(2*i, 1) = -rocket.lbx(x_indexes(i));
            end
            
            [K,Qf,~] = dlqr(mpc.A,mpc.B,Q,R);
            K = -K; 
            %P = dlyap(mpc.A,Q);
            P = Qf;
            %% Set up the MPC cost and constraints using the computed set-point            
            con = (X(:,2) == mpc.A*X(:,1) + mpc.B*U(:,1)) + (M*U(:,1) <= m);
            obj = (U(:,1) - u_ref) *R* (U(:,1) - u_ref);
            for i = 2:N-1
                con = con + (X(:,i+1) == mpc.A*X(:,i) + mpc.B*U(:,i));
                con = con + (M*U(:,i) <= m) + (F*X(:,i) <= f); 
                obj = obj + (X(:,i)-x_ref)'*Q*(X(:,i)-x_ref) + ...
                            (U(:,i)-u_ref)'*R*(U(:,i)-u_ref);
            end
            obj = obj + (X(:,N)-x_ref)'*P*(X(:,N)-x_ref);
            
            %% Optimize
            diagnostics = solvesdp(con, obj, sdpsettings('verbose',0));
            
            if diagnostics.problem == 0
            % Good! 
            else
                throw(MException('',yalmiperror(diagnostics.problem)));
            end
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Return YALMIP optimizer object
            ctrl_opti = optimizer(con, obj, sdpsettings('solver','gurobi'), ...
                {X(:,1), x_ref, u_ref}, U(:,1));
        end
        
        % Design a YALMIP optimizer object that takes a position reference
        % and returns a feasible steady-state state and input (xs, us)
        function target_opti = setup_steady_state_target(mpc)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % INPUTS
            %   ref    - reference to track
            % OUTPUTS
            %   xs, us - steady-state target
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            nx = size(mpc.A, 1);

            % Steady-state targets
            xs = sdpvar(nx, 1);
            us = sdpvar;
            
            % Reference position (Ignore this before Todo 3.2)
            ref = sdpvar;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            % You can use the matrices mpc.A, mpc.B, mpc.C and mpc.D            
            global rocket
            u_index = 2;            
            x_indexes = [2,5,7,10];
            
            % u in U = { u | Mu <= m }
            M = [1;-1]; m = [rocket.ubu(u_index); -rocket.lbu(u_index)];
            % x in X = { x | Fx <= f }
            F = zeros(nx*2, nx);
            f = zeros(nx*2, 1);
            for i = 1:nx
                F(2*i - 1, i) = 1;
                F(2*i , i) = -1;
                f(2*i-1, 1) = rocket.ubx(x_indexes(i));
                f(2*i, 1) = -rocket.lbx(x_indexes(i));
            end
            
            con = [M*us <= m, F*xs <= f, ...
            xs == mpc.A*xs + mpc.B*us    ,...
            ref == mpc.C*xs ];

            obj = us^2;
            diagnostics = solvesdp(con,obj,sdpsettings('verbose',0));

            if diagnostics.problem == 0
            % Good! 
            elseif diagnostics.problem == 1
            throw(MException('','Infeasible'));
            % TODO should we find a closest solution?
            else
            throw(MException('','Something else happened'));
            end
            
            % YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE YOUR CODE HERE
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Compute the steady-state target
            target_opti = optimizer(con, obj, sdpsettings('solver', 'gurobi'), ref, {xs, us});
        end
    end
end
