%% compute_trim.m
clear; clc;

aerosonde_parameters;

gamma = 0*pi/180;
R     = 10000000000;
Va    = 200;

%% --- Build 12-state Euler initial condition ---
x0_12 = [0; 0; -200; Va; 0; 0; 0; 0; 0; 0; 0; 0];

%% --- Convert to 13-state quaternion for Simulink ---
q0 = Euler2Quaternion(x0_12(7), x0_12(8), x0_12(9));    % ← 3 separate args
x0_13 = [x0_12(1:6); q0; x0_12(10:12)];

%% --- Inputs / outputs ---
u0 = [0; 0; 0; 1];
iu = [];
ix = [];
y0 = [Va; 0; 0];
iy = [1, 3];

%% --- Derivatives (13-state) ---
dx0_13 = [0; 0; -Va*sin(gamma); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
if R ~= Inf, dx0_13(13) = Va*cos(gamma)/R; end

idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13];

%% --- Run trim ---
[x_trim13, u_trim, y_trim, dx_trim13] = ...
    trim('mavsim_trim', x0_13, u0, y0, ix, iu, iy, dx0_13, idx);

%% --- Convert result back to 12-state Euler ---
[phi, theta, psi] = Quaternion2Euler(x_trim13(7:10));
x_trim = [x_trim13(1:6); phi; theta; psi; x_trim13(11:13)];

%% --- Sanity check ---
dx_trim12 = [dx_trim13(1:6); zeros(3,1); dx_trim13(11:13)];
dx0_12    = [dx0_13(1:6);    zeros(3,1); dx0_13(11:13)];
fprintf('Trim residual: %.6e\n', norm(dx_trim12(3:end) - dx0_12(3:end)));

%% --- Set MAV initial conditions ---
MAV.pn0    = 0;
MAV.pe0    = 0;
MAV.pd0    = -200;
MAV.u0     = x_trim(4);
MAV.v0     = x_trim(5);
MAV.w0     = x_trim(6);
MAV.phi0   = x_trim(7);
MAV.theta0 = x_trim(8);
MAV.psi0   = x_trim(9);
MAV.p0     = x_trim(10);
MAV.q0     = x_trim(11);
MAV.r0     = x_trim(12);

fprintf('Trim complete:\n');
fprintf('  alpha   = %.4f deg\n', y_trim(2)*180/pi);
fprintf('  theta   = %.4f deg\n', x_trim(8)*180/pi);
fprintf('  delta_e = %.4f  delta_t = %.4f\n', u_trim(1), u_trim(4));
