addpath(fullfile('..', 'src'));
rmpath('Deliverable_3_1\');
rmpath('Deliverable_3_2\');
rmpath('Deliverable_4_1\');
addpath('Deliverable_5_1/');
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
%% 
mpc_y = MPC_Control_y(sys_y, Ts, H);
%% 
mpc_z = MPC_Control_z(sys_z, Ts, H);
%% 
mpc_roll = MPC_Control_roll(sys_roll, Ts, H);
%% Merge
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% Setup reference function
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);
%%
% Tweak rocket mass
rocket.mass = 1.783;%2.0; % 1.8;

[T, X, U, Ref] = rocket.simulate_f(x0, Tf, mpc, ref);
%%
% Plot pose
rocket.anim_rate = 10;
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation with mass disturbance impact';
%%
mpc_z = MPC_Control_z(sys_z, Ts, H);
mpc = rocket.merge_lin_controllers(xs, us, mpc_x, mpc_y, mpc_z, mpc_roll);
% Setup reference function
Tf = 30;
ref = @(t_, x_) rocket.MPC_ref(t_, Tf);
x0 = zeros(12,1);
% Tweak rocket mass
rocket.mass = 1.783;

[T, X, U, Ref, Z_hat] = rocket.simulate_f_est_z(x0, Tf, mpc, ref, mpc_z, sys_z);

%%
% Plot pose
rocket.anim_rate = 5;
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation with mass disturbance impact with offset-free';

%%
% Tweak rocket mass
rocket.mass = 1.8;

[T, X, U, Ref, Z_hat] = rocket.simulate_f_est_z(x0, Tf, mpc, ref, mpc_z, sys_z);

rocket.anim_rate = 5;
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation with mass disturbance impact with offset-free';
%%
% Tweak rocket mass
rocket.mass = 2.0;

[T, X, U, Ref, Z_hat] = rocket.simulate_f_est_z(x0, Tf, mpc, ref, mpc_z, sys_z);

rocket.anim_rate = 5;
ph = rocket.plotvis(T, X, U, Ref);
ph.fig.Name = 'Merged lin. MPC in nonlinear simulation with mass disturbance impact with offset-free';

