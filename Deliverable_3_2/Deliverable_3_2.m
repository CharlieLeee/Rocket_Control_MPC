addpath(fullfile('..', 'src'));
rmpath('Deliverable_3_1\');
addpath('Deliverable_3_2\');
clc; clear;
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
x0 = zeros(12,1); %x0([10,11,12]) = 5; x0(6) = pi/4;

%% 
mpc_x = MPC_Control_x(sys_x, Ts, H);
x_ref = -5;
[T, X_sub, U_sub] = rocket.simulate(sys_x, x0(x_indexes_list{1}), Tf, @mpc_x.get_u, x_ref);
ph_x = rocket.plotvis_sub(T, X_sub, U_sub, sys_x, xs, us, x_ref);
%% 
mpc_y = MPC_Control_y(sys_y, Ts, H);
y_ref = -5;
[T, Y_sub, U_sub] = rocket.simulate(sys_y, x0(x_indexes_list{2}), Tf, @mpc_y.get_u, y_ref);
ph_y = rocket.plotvis_sub(T, Y_sub, U_sub, sys_y, xs, us, y_ref);
%% 
mpc_z = MPC_Control_z(sys_z, Ts, H);
z_ref = -5;
[T, Z_sub, U_sub] = rocket.simulate(sys_z, x0(x_indexes_list{3}), Tf, @mpc_z.get_u, z_ref);
ph_z = rocket.plotvis_sub(T, Z_sub, U_sub, sys_z, xs, us, z_ref);
%% 
mpc_r = MPC_Control_roll(sys_roll, Ts, H);
r_ref = pi/4;  
[T, R_sub, U_sub] = rocket.simulate(sys_roll, x0(x_indexes_list{4}), Tf, @mpc_r.get_u, r_ref);
ph_r = rocket.plotvis_sub(T, R_sub, U_sub, sys_roll, xs, us);