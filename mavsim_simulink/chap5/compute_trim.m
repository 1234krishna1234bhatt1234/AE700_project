% compute_trim.m
% Computes trim (equilibrium) conditions for the jet aircraft model.
%
% Inputs (edit these before running):
%   Va    - desired airspeed [m/s]
%   gamma - desired flight-path angle [rad]  (+up, -down)
%   R     - desired turn radius [m]  (+right, -left, Inf = straight)
%
% Prerequisites: aerosonde_parameters.m must have been run so that P and MAV
%                exist in the workspace.
%
% The trim model 'mavsim_trim' must be a 12-state Euler-angle Simulink model
% with:
%   Inputs : [delta_e; delta_a; delta_r; delta_t]
%   Outputs: [Va; alpha; beta]   (first three outputs used by iy)
%
% -------------------------------------------------------------------------

gamma = 0 * pi/180;    % desired flight-path angle (rad)
R     = Inf;           % turn radius (m): Inf = straight, +R right, -R left
Va    = 210;           % desired airspeed (m/s)  ~203 m/s gives alpha~6 deg

% --- Initial state guess ---
alpha_guess = 6 * pi/180;        % AoA initial guess [rad]
theta0      = alpha_guess + gamma; % pitch = AoA + gamma for coordinated flight
phi0        = 0;                  % wings level

x0 = [0;                     % 1  pn [m]
      0;                     % 2  pe [m]
      -10000;                % 3  pd [m]  (10 km cruise altitude)
      Va*cos(alpha_guess);   % 4  u [m/s]
      0;                     % 5  v [m/s]
      Va*sin(alpha_guess);   % 6  w [m/s]
      phi0;                  % 7  phi [rad]
      theta0;                % 8  theta [rad]
      0;                     % 9  psi [rad]
      0;                     % 10 p [rad/s]
      0;                     % 11 q [rad/s]
      0];                    % 12 r [rad/s]

% specify which states to hold equal to x0 (none -- let trim find all)
ix = [];

% --- Initial control guess ---
qbar0     = 0.5 * P.rho * Va^2;
CD_guess  = P.C_D_0 + P.C_D_alpha * alpha_guess;
Drag_est  = qbar0 * P.S_wing * CD_guess;
dt0       = min(1, max(0.01, Drag_est / P.T_max));  % throttle fraction estimate

u0 = [0;    % 1 delta_e [rad]
      0;    % 2 delta_a [rad]
      0;    % 3 delta_r [rad]
      dt0]; % 4 delta_t [0-1]

% specify which inputs to hold constant (none)
iu = [];

% --- Desired outputs ---
y0 = [Va;   % 1 Va [m/s]   -- hold airspeed
      0;    % 2 alpha      -- not constrained
      0];   % 3 beta [rad] -- hold sideslip at zero

% hold Va (output 1) and beta (output 3) constant
iy = [1, 3];

% --- Desired state derivatives (trim = equilibrium) ---
dx0 = zeros(12, 1);
% for a turn: non-zero psi_dot
if ~isinf(R) && R ~= 0
    dx0(9) = Va * cos(gamma) / R;  % psi_dot [rad/s]
end

% hold all derivatives from pd_dot (3) through r_dot (12) at dx0 values
idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12];

% --- Run trim algorithm ---
[x_trim, u_trim, y_trim, dx_trim] = trim('mavsim_trim', x0, u0, y0, ix, iu, iy, dx0, idx);

% Check residual (should be very small, < 1e-6)
fprintf('Trim residual norm: %.6e  (should be < 1e-6)\n', norm(dx_trim(3:end) - dx0(3:end)));

% --- Extract key trim quantities for TF / SS scripts ---
Va_trim    = y_trim(1);           % [m/s]
alpha_trim = y_trim(2);           % [rad]
beta_trim  = y_trim(3);           % [rad]  (should be ~0)
theta_trim = x_trim(8);           % [rad]
phi_trim   = x_trim(7);           % [rad]
gamma_trim = gamma;               % [rad]

fprintf('Va_trim    = %.4f m/s\n',  Va_trim);
fprintf('alpha_trim = %.4f rad  (%.2f deg)\n', alpha_trim, alpha_trim*180/pi);
fprintf('theta_trim = %.4f rad  (%.2f deg)\n', theta_trim, theta_trim*180/pi);
fprintf('delta_e_trim = %.6f rad\n', u_trim(1));
fprintf('delta_t_trim = %.6f\n',    u_trim(4));

% --- Set initial conditions for main quaternion simulation ---
MAV.pn0    = 0;
MAV.pe0    = 0;
MAV.pd0    = -10000;             % altitude = 10 km
MAV.u0     = x_trim(4);         % body-x velocity
MAV.v0     = x_trim(5);         % body-y velocity
MAV.w0     = x_trim(6);         % body-z velocity
MAV.phi0   = x_trim(7);         % roll
MAV.theta0 = x_trim(8);         % pitch
MAV.psi0   = x_trim(9);         % heading

% Convert Euler trim angles to quaternion for quaternion-based main simulation
e_q = Euler2Quaternion(x_trim(7), x_trim(8), x_trim(9));
MAV.e0 = e_q(1);
MAV.e1 = e_q(2);
MAV.e2 = e_q(3);
MAV.e3 = e_q(4);

MAV.p0  = x_trim(10);           % roll rate
MAV.q0  = x_trim(11);           % pitch rate
MAV.r0  = x_trim(12);           % yaw rate

% Store trim inputs for autopilot / TF / SS scripts
P.u_trim = u_trim;
P.x_trim = x_trim;
