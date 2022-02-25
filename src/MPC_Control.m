classdef MPC_Control
    properties
        ctrl_opti   % YALMIP optimizer object to compute control law
        target_opti % YALMIP optimizer object to compute steady-state target
        
        % Discrete-time system matrices
        A, B, C, D
    end
    
    methods
        function mpc = MPC_Control(sys, Ts, H)
            
            % Discretize the system and extract the A,B,C,D matrices
            sys_d = c2d(sys, Ts);
            [mpc.A, mpc.B, mpc.C, mpc.D] = ssdata(sys_d);
            
            % Horizon
            mpc.ctrl_opti   = mpc.setup_controller(Ts, H);
            mpc.target_opti = mpc.setup_steady_state_target();
        end
        
        % Compute the MPC controller
        function u = get_u(mpc, x, ref)
            
            % Compute the target ...
            if nargin >= 3
                % ... from reference
                if length(struct(mpc.target_opti).diminOrig) == 2
                    if length(x) == 3
                        d_est = x(end);
                    else
                        d_est = 0;
                    end
                    [target, solve_status] = mpc.target_opti(ref, d_est);
                else
                    [target, solve_status] = mpc.target_opti(ref);
                end
                [ref_x, ref_u] = deal(target{:});
                if solve_status ~= 0
                    solve_status_str = yalmiperror(solve_status);
                    fprintf([' [Target: ' solve_status_str(1:end-1) '] ']);
                    u = nan(struct(mpc.ctrl_opti).dimoutOrig{1}); return
                end
            else
                % ... set origin
                ref_x = zeros(size(mpc.A,1),1);
                ref_u = 0;
            end
            
            % Compute the control action
            if length(struct(mpc.ctrl_opti).diminOrig) == 4
                if length(x) == 3
                    d_est = x(end);
                    x = x(1:end-1);
                else
                    d_est = 0;
                end
                [u, solve_status] = mpc.ctrl_opti({x, ref_x, ref_u, d_est});
            else
                [u, solve_status] = mpc.ctrl_opti({x, ref_x, ref_u});
            end
            if solve_status ~= 0
                solve_status_str = yalmiperror(solve_status);
                fprintf([' [Control: ' solve_status_str(1:end-1) '] ']);
                u = nan(size(u));
            end
            
        end
    end
    
    methods (Abstract)
        setup_controller(mpc)
        setup_steady_state_target(mpc)
    end
end
