%% 1. Setup Parameters
aerosonde_parameters; 

% Flight Condition
gamma = 0 * pi/180;    % Level flight
R     = Inf;           % Straight line
Va    = 200;           % 100 m/s for the 747

%% 2. Set Initial Conditions (12 States)
% [pn; pe; pd; u; v; w; phi; theta; psi; p; q; r]
x0 = [0; 0; -500; Va; 0; 0; 0; gamma; 0; 0; 0; 0];
ix = []; % No states are fixed to x0 exactly

%% 3. Set Initial Inputs
% [delta_e; delta_a; delta_r; delta_t]
u0 = [0; 0; 0; 0.5]; % Start with 50% throttle
iu = [];

%% 4. Define Constant Outputs
% This must match the size of the Outport in your Simulink model!
% If your model outputs [Va; alpha; beta], then y0 is length 3.
y0 = [Va; gamma; 0]; 
iy = [1, 3]; % Constrain Va (1) and Side-slip (3)

%% 5. Define Constant Derivatives
% Standard steady-state: derivatives = 0, except for altitude change if climbing.
dx0 = [0; 0; -Va*sin(gamma); 0; 0; 0; 0; 0; 0; 0; 0; 0];

if R ~= Inf
    dx0(9) = Va*cos(gamma)/R; % Target psidot for a turn
end

% idx: Derivatives to hold constant. 
% We lock 3-12 (Velocity rates and Angular rates). 
% We leave 1 & 2 (North/East position) free.
idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12];

%% 6. Compute Trim
% Run the solver using the 9-argument format
[x_trim, u_trim, y_trim, dx_trim] = trim('mavsim_trim', x0, u0, y0, ix, iu, iy, dx0, idx);

% Success Check
residual = norm(dx_trim(idx) - dx0(idx));
fprintf('Trim Residual: %.2e\n', residual);

%% 7. Update MAV Initial Conditions
MAV.pn0    = 0;
MAV.pe0    = 0;
MAV.pd0    = -500; % Starting Altitude
MAV.u0     = x_trim(4);
MAV.v0     = x_trim(5);
MAV.w0     = x_trim(6);
MAV.phi0   = x_trim(7);
MAV.theta0 = x_trim(8);
MAV.psi0   = x_trim(9);
MAV.p0     = x_trim(10);
MAV.q0     = x_trim(11);
MAV.r0     = x_trim(12);

fprintf('Trimmed Inputs: [Elev: %.2f deg] | [Ail: %.2f deg] | [Rud: %.2f deg] | [Thr: %.2f%%]\n', u_trim(1)*180/pi, u_trim(2)*180/pi, u_trim(3)*180/pi, u_trim(4)*100);
% %% compute_trim.m
% 
% aerosonde_parameters;
% 
% gamma = 0*pi/180;
% R     = 10000000;
% Va    = 200;
% 
% %% --- Build 12-state Euler initial condition ---
% x0_12 = [0; 0; -200; Va; 0; 0; 0; 2*pi/180; 0; 0; 0; 0];
% 
% %% --- Convert to 13-state quaternion for Simulink ---
% q0 = Euler2Quaternion(x0_12(7), x0_12(8), x0_12(9));    % ← 3 separate args
% x0_13 = [x0_12(1:6); q0; x0_12(10:12)];
% 
% %% --- Inputs / outputs ---
% u0 = [0; 0; 0; 0.5];
% iu = [];
% ix = [];
% y0 = [Va; 0; 0];
% iy = [1, 3];
% 
% %% --- Derivatives (13-state) ---
% dx0_13 = [0; 0; -Va*sin(gamma); 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% if R ~= Inf, dx0_13(13) = Va*cos(gamma)/R; end
% 
% idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13];
% 
% %% --- Run trim ---
% [x_trim13, u_trim, y_trim, dx_trim13] = ...
%     trim('mavsim_trim', x0_13, u0, y0, ix, iu, iy, dx0_13, idx);
% 
% %% --- Convert result back to 12-state Euler ---
% [phi, theta, psi] = Quaternion2Euler(x_trim13(7:10));
% x_trim = [x_trim13(1:6); phi; theta; psi; x_trim13(11:13)];
% 
% %% --- Sanity check ---
% dx_trim12 = [dx_trim13(1:6); zeros(3,1); dx_trim13(11:13)];
% dx0_12    = [dx0_13(1:6);    zeros(3,1); dx0_13(11:13)];
% fprintf('Trim residual: %.6e\n', norm(dx_trim12(3:end) - dx0_12(3:end)));
% 
% %% --- Set MAV initial conditions as trim conditions ---
% MAV.pn0    = x_trim(1);
% MAV.pe0    = x_trim(2);
% MAV.pd0    = x_trim(3);
% MAV.u0     = x_trim(4);
% MAV.v0     = x_trim(5);
% MAV.w0     = x_trim(6);
% MAV.phi0   = x_trim(7);
% MAV.theta0 = x_trim(8);
% MAV.psi0   = x_trim(9);
% MAV.p0     = x_trim(10);
% MAV.q0     = x_trim(11);
% MAV.r0     = x_trim(12);
% 
% fprintf('Trim complete:\n');
% fprintf('  alpha   = %.4f deg\n', y_trim(2)*180/pi);
% fprintf('  theta   = %.4f deg\n', x_trim(8)*180/pi);
% fprintf('  delta_e = %.4f delta_a = %.4f delta_r = %.4f delta_t = %.4f\n', u_trim(1), u_trim(2),u_trim(3),u_trim(4));



