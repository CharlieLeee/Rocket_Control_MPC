addpath(fullfile('..', 'src'));

%% TODO: This file should produce all the plots for the deliverable
Ts = 1/20; % Sample time 
global rocket
rocket = Rocket(Ts);
[xs, us] = rocket.trim(); 
sys = rocket.linearize(xs, us); 
[sys_x, sys_y, sys_z, sys_roll] = rocket.decompose(sys, xs, us);

H = 5; % Horizon length in seconds 
x_indexes_list = {[2,5,7,10],
                  [1,4,8,11],
                  [9,12],
                  [3,6]};                  
Tf = 10;
x0 = zeros(12,1); x0([10,11,12]) = 5; x0(6) = pi/4;

%% X
mpc_x = MPC_Control_x(sys_x, Ts, H);
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0(x_indexes_list{1}), Tf, @mpc_x.get_u, 0);
ph_x = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us);
%% Y
mpc_y = MPC_Control_y(sys_y, Ts, H);
[T, Y_sub, U_sub] = rocket.simulate(sys_y, x0(x_indexes_list{2}), Tf, @mpc_y.get_u, 0);
ph_y = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us);
%% Z
mpc_z = MPC_Control_z(sys_z, Ts, H);
[T, Z_sub, U_sub] = rocket.simulate(sys_z, x0(x_indexes_list{3}), Tf, @mpc_z.get_u, 0);
ph_z = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us);
%% Roll
mpc_r = MPC_Control_roll(sys_roll, Ts, H);
[T, R_sub, U_sub] = rocket.simulate(sys_roll, x0(x_indexes_list{4}), Tf, @mpc_r.get_u, 0);
ph_r = rocket.plotvis_sub(T, R_sub, U_sub, sys_roll, xs, us);