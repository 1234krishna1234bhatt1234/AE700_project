% compute_tf_model.m
% x_trim, u_trim, y_trim must be in workspace from compute_trim.m
% P must be in workspace (aircraft parameters, same as forces_moments.m)
% If you used MAV struct: P = MAV;

%% --- Step 1: Extract trim values ---
Va_trim    = y_trim(1);          % trimmed airspeed (m/s)
alpha_trim = y_trim(2);          % trimmed angle of attack (rad)
beta_trim  = y_trim(3);          % trimmed sideslip (rad), should be ~0
theta_trim = x_trim(8);          % trimmed pitch angle (rad)
phi_trim   = x_trim(7);          % trimmed roll angle (rad)
delta_e_trim = u_trim(1);        % trimmed elevator
delta_t_trim = u_trim(4);        % trimmed throttle

%% --- Step 2: Compute transfer function coefficients ---
% All formulas from Beard & McLain, "Small Unmanned Aircraft", Ch. 5

% ---- ROLL subsystem (B&M Eq. 5.38) ----
% phi_ddot = a_phi1*phi_dot + a_phi2*delta_a
% TF: phi(s)/delta_a(s) = a_phi2 / (s^2 + a_phi1*s)
%
% a_phi1 = -[Γ3·C_l_p + Γ4·C_n_p] · (ρ·Va·S·b²/4)
% a_phi2 = [Γ3·C_l_δa + Γ4·C_n_δa] · (ρ·Va²·S·b/2)
a_phi1 = -(P.Gamma3 * P.C_ell_p + P.Gamma4 * P.C_n_p) ...
          * (P.rho * Va_trim * P.S_wing * P.b^2 / 4);
a_phi2 =  (P.Gamma3 * P.C_ell_delta_a + P.Gamma4 * P.C_n_delta_a) ...
          * (P.rho * Va_trim^2 * P.S_wing * P.b / 2);

% ---- PITCH subsystem (B&M Eq. 5.23–5.25) ----
% theta_ddot = a_theta1*theta_dot + a_theta2*theta + a_theta3*delta_e
% TF: theta(s)/delta_e(s) = a_theta3 / (s^2 + a_theta1*s + a_theta2)
%
% a_theta1 = -(ρ·Va·S·c²)/(4·Jy) · C_m_q         [pitch damping, >0 since C_m_q<0]
% a_theta2 = -(ρ·Va²·S·c)/(2·Jy) · C_m_alpha      [pitch stiffness, >0 since C_m_α<0]
% a_theta3 = -(ρ·Va²·S·c)/(2·Jy) · C_m_delta_e    [elevator effectiveness]
a_theta1 = -(P.rho * Va_trim * P.S_wing * P.c^2) / (4 * P.Jy) * P.C_m_q;
a_theta2 = -(P.rho * Va_trim^2 * P.S_wing * P.c) / (2 * P.Jy) * P.C_m_alpha;
a_theta3 = -(P.rho * Va_trim^2 * P.S_wing * P.c) / (2 * P.Jy) * P.C_m_delta_e;

% DC gain of the closed pitch inner loop (used in altitude loop design):
%   K_theta_DC = a_theta3*kp_pitch / (a_theta2 + a_theta3*kp_pitch)
% (computed in compute_autopilot_gains.m)

% ---- AIRSPEED subsystem (B&M Eq. 5.13–5.15) ----
% Va_dot = -a_V1*(Va - Va*) + a_V2*(delta_t - delta_t*) - a_V3*(theta - theta*)
% TF throttle: Va(s)/delta_t(s) = a_V2 / (s + a_V1)
% TF pitch:    Va(s)/theta(s)   = -a_V3 / (s + a_V1)
%
% Propulsion model (from forces_moments.m):
%   Thrust = 0.5·ρ·S_prop·C_prop·[(k_motor·δt)² - Va²]
%   → ∂T/∂Va     = -ρ·S_prop·C_prop·Va_trim
%   → ∂T/∂delta_t = ρ·S_prop·C_prop·k_motor²·delta_t_trim
%
% Aerodynamic drag gradient:
%   C_D_trim = C_D_0 + C_D_alpha*alpha_trim + C_D_delta_e*delta_e_trim
%   ∂D/∂Va = ρ·Va·S·C_D_trim
%
% a_V1 = (1/m)·[∂D/∂Va - ∂T/∂Va]
%       = (ρ·Va_trim/m)·[S·C_D_trim + S_prop·C_prop]
C_D_trim = P.C_D_0 + P.C_D_alpha * alpha_trim + P.C_D_delta_e * delta_e_trim;
% a_V1 = (P.rho * Va_trim / P.mass) * (P.S_wing * C_D_trim + P.S_prop * P.C_prop);
a_V1 =0;
% a_V2 = (1/m)·∂T/∂delta_t = ρ·S_prop·C_prop·k_motor²·delta_t_trim / m
a_V2 = P.T_max / P.mass;

% a_V3 = g·cos(gamma_trim) ≈ g   (gravity coupling through pitch angle)
a_V3 = P.gravity;

% ---- SIDESLIP subsystem (B&M Eq. 5.36) ----
% beta_dot = a_beta1*beta + a_beta2*delta_r
% TF: beta(s)/delta_r(s) = a_beta2 / (s + a_beta1)
%
% a_beta1 = -(ρ·Va·S)/(2m) · C_Y_beta    [>0 since C_Y_beta<0]
% a_beta2 =  (ρ·Va·S)/(2m) · C_Y_delta_r [sign depends on convention]
a_beta1 = -(P.rho * Va_trim * P.S_wing) / (2 * P.mass) * P.C_Y_beta;
a_beta2 =  (P.rho * Va_trim * P.S_wing) / (2 * P.mass) * P.C_Y_delta_r;

%% --- Step 3: Define transfer functions ---
T_phi_delta_a   = tf([a_phi2],         [1, a_phi1, 0]);
T_chi_phi       = tf([P.gravity/Va_trim], [1, 0]);
T_theta_delta_e = tf(a_theta3,          [1, a_theta1, a_theta2]);
T_h_theta       = tf([Va_trim],         [1, 0]);
T_h_Va          = tf([theta_trim],      [1, 0]);
T_Va_delta_t    = tf([a_V2],            [1, a_V1]);
T_Va_theta      = tf([-a_V3],           [1, a_V1]);
T_v_delta_r     = tf([a_beta2],         [1, a_beta1]);

%% --- Step 4: Save coefficients for autopilot design ---
Va = Va_trim;  % alias for compute_autopilot_gains.m
save('transfer_function_coef', ...
    'a_phi1','a_phi2', ...
    'a_theta1','a_theta2','a_theta3', ...
    'a_V1','a_V2','a_V3', ...
    'a_beta1','a_beta2', ...
    'Va_trim','Va','theta_trim','alpha_trim');

disp('Transfer function coefficients:')
disp(['  a_phi1 = ', num2str(a_phi1), '  a_phi2 = ', num2str(a_phi2)])
disp(['  a_theta1 = ', num2str(a_theta1), '  a_theta2 = ', num2str(a_theta2), '  a_theta3 = ', num2str(a_theta3)])
disp(['  a_V1 = ', num2str(a_V1), '  a_V2 = ', num2str(a_V2), '  a_V3 = ', num2str(a_V3)])
disp(['  a_beta1 = ', num2str(a_beta1), '  a_beta2 = ', num2str(a_beta2)])
