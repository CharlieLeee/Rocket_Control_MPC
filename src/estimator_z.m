function z_hat_next = estimator_z(z_hat, u, x, mpc_z, sys_z)
%Z_ESTIMATOR Estimation function for full state vector + disturbance
%   This estimation function wraps the state and disturbance estimation
%   into a more general estimation function for the full state augmented by
%   additional estimation variables.
%   z_hat_next = [x_hat; est] contains a full state estimate x_hat and
%   additional estimated variables est.
%   The entries of x_hat that belong to the z sub-system contain the state
%   estimate of the z sub-system observer. The remaining entries contain
%   the true states.

% Indices of the z sub-system
idx_z = sys_z.UserData.idx;
idu_z = sys_z.UserData.idu;
idy_z = sys_z.UserData.idy;
us_z  = sys_z.UserData.us;

% Estimation equation
z_hat_z = mpc_z.A_bar * [z_hat(idx_z); z_hat(13:end)] + ...
    mpc_z.B_bar * (u(idu_z) - us_z) + ...
    mpc_z.L * (mpc_z.C_bar * [z_hat(idx_z); z_hat(13:end)] - x(idy_z));

% Hand through true states except for estimated ones
x_hat = x;
x_hat(idx_z) = z_hat_z(1:2);

% Concatenate estimation variable z = [x; est]
z_hat_next = [x_hat; z_hat_z(3)];

end