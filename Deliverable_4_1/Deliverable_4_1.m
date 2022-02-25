addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable
Ts = 0.05;
rocket = Rocket(Ts); 
[xs, us] = rocket.trim(); 
sys = rocket.linearize(xs, us);
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

H = 5;
mpc_x = MPC_Control_x(sys_x, Ts, H);
mpc_y = MPC_Control_y(sys_y, Ts, H);
mpc_z = MPC_Control_z(sys_z, Ts, H);
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);

% Merge four sub−system controllers into one full−system controller 
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% Setup reference function 
Tf = 30; 
ref = @(t_ , x_ ) rocket.MPC_ref(t_ , Tf);

x0 = zeros(12,1); 
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);
% Plot pose 
rocket.anim_rate = 5; % Increase this to make the animation faster 
ph = rocket.plotvis(T, X, U, Ref); 
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation'; % Set a figure title