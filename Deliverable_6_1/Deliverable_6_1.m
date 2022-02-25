addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/10; % Note that we choose a larger Ts here to speed up the simulation 
rocket = Rocket(Ts);
H = 1;
nmpc = NMPC_Control(rocket, H);
x0 = zeros(12,1); 

%%
% MPC reference with default maximum roll = 15 deg 
Tf = 30; 
ref = @(t , x ) rocket.MPC_ref(t , Tf);

%%
% MPC reference with specified maximum roll = 50 deg 
Tf = 30; 
roll_max = deg2rad(50); 
ref = @(t , x ) rocket.MPC_ref(t , Tf, roll_max);

%%
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, nmpc, ref);
%%
rocket.anim_rate = 5.0;
ph = rocket.plotvis(T, X, U, Ref);

%%
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

lmpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
%%
[T, X, U, Ref] = rocket.simulate_f(x0, Tf, lmpc, ref);