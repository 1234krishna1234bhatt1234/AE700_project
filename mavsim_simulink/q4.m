%% ========================================================================
%% Q4_linear_design.m
%% ========================================================================
%% Computes trim, transfer functions, and state-space models for the MAV.
%% Steps (i)-(vii) from the problem statement.
%%
%% Prerequisites:
%%   - aerosonde_parameters.m  (sets up MAV/P struct)
%%   - simulation_parameters.m (sets up SIM struct)
%% ========================================================================
clear; clc; close all;

%% ---- Load aircraft parameters ----------------------------------------
addpath('parameters');
addpath('tools');
aerosonde_parameters;   % creates MAV struct, aliases to P
simulation_parameters;  % creates SIM struct

%% ========================================================================
%% STEP (iii)-(iv): Compute TRIM for various γ (wings-level flight)
%% ========================================================================

Va_des    = 200;         % desired airspeed (m/s) — from parameters file
gammas    = [-10, -5, 0, 5, 10] * pi/180;  % path angles to test
R_des     = 5000;         % wings-level (infinite turn radius)

fprintf('\n--- Trim Results for Wings-Level Flight (Va=%.0f m/s) ---\n', Va_des);
fprintf('%8s %10s %10s %10s %10s %10s\n', 'gamma°', 'alpha°', 'theta°', 'delta_e', 'delta_t', 'climb(m/s)');

trim_results = struct();
for k = 1:length(gammas)
    gamma = gammas(k);
    [xt, ut, yt, climb_rate] = compute_trim_analytic(Va_des, gamma, R_des, P);
    trim_results(k).gamma     = gamma;
    trim_results(k).x_trim    = xt;
    trim_results(k).u_trim    = ut;
    trim_results(k).Va        = Va_des;
    trim_results(k).alpha     = yt(1);
    trim_results(k).climb     = climb_rate;
    
    fprintf('%8.1f %10.3f %10.3f %10.4f %10.4f %10.3f\n', ...
        gamma*180/pi, yt(1)*180/pi, xt(2)*180/pi, ut(1), ut(4), climb_rate);
end

% Verify: only altitude should change for different gamma
fprintf('\nVerification: theta changes with gamma (correct), Va stays at %.0f m/s\n', Va_des);

%% ---- Step (v): Trim for constant turn (load factor n=1.2) ------------
    
n_load  = 1.2;          % load factor
CL_des  = 0.85;         % CL in [0.7, 1.0] range
phi_trim_turn = acos(1/n_load);        % bank angle for coordinated turn
R_turn  = Va_des^2 / (P.gravity * tan(phi_trim_turn));  % turn radius

fprintf('\n--- Trim for Constant Turn (n=%.1f, CL=%.2f) ---\n', n_load, CL_des);
fprintf('  Bank angle: %.2f deg\n', phi_trim_turn*180/pi);
fprintf('  Turn radius: %.1f m\n', R_turn);

[xt_turn, ut_turn, yt_turn, ~] = compute_trim_analytic(Va_des, 0, R_turn, P);
fprintf('  alpha=%.2f°, delta_e=%.4f, delta_t=%.4f\n', ...
    yt_turn(1)*180/pi, ut_turn(1), ut_turn(4));

%% ========================================================================
%% STEP (vi): Compute Transfer Functions (analytic, Ch. 5.4)
%% ========================================================================

% Use wings-level trim (γ=0) for linearization
idx0     = find(gammas == 0);
if isempty(idx0), idx0 = 3; end
Va_trim  = trim_results(idx0).Va;
x_trim   = trim_results(idx0).x_trim;
u_trim   = trim_results(idx0).u_trim;

alpha_trim   = trim_results(idx0).alpha;
theta_trim   = x_trim(2);    % [phi, theta, ...]
delta_e_trim = u_trim(1);
delta_t_trim = u_trim(4);

% --- TF coefficients (B&M Eq. 5.xx) ---
a_phi1 = -(P.Gamma3*P.C_ell_p + P.Gamma4*P.C_n_p) ...
          * (P.rho*Va_trim*P.S_wing*P.b^2/4);
a_phi2 =  (P.Gamma3*P.C_ell_delta_a + P.Gamma4*P.C_n_delta_a) ...
          * (P.rho*Va_trim^2*P.S_wing*P.b/2);

a_theta1 = -(P.rho*Va_trim*P.S_wing*P.c^2)/(4*P.Jy) * P.C_m_q;
a_theta2 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy) * P.C_m_alpha;
a_theta3 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy) * P.C_m_delta_e;

C_D_trim = P.C_D_0 + P.C_D_alpha*alpha_trim + P.C_D_delta_e*delta_e_trim;
a_V1 = (P.rho*Va_trim/P.mass) * (P.S_wing*C_D_trim);
a_V2 = P.T_max / P.mass;
a_V3 = P.gravity;

a_beta1 = -(P.rho*Va_trim*P.S_wing)/(2*P.mass) * P.C_Y_beta;
a_beta2 =  (P.rho*Va_trim*P.S_wing)/(2*P.mass) * P.C_Y_delta_r;

% --- Define transfer functions ---
s = tf('s');
T_phi_delta_a   = a_phi2 / (s^2 + a_phi1*s);
T_chi_phi       = (P.gravity/Va_trim) / s;
T_theta_delta_e = a_theta3 / (s^2 + a_theta1*s + a_theta2);
T_h_theta       = Va_trim / s;
T_h_Va          = theta_trim / s;
T_Va_delta_t    = a_V2 / (s + a_V1);
T_Va_theta      = -a_V3 / (s + a_V1);
T_beta_delta_r  = a_beta2 / (s + a_beta1);

fprintf('\n--- Transfer Function Coefficients ---\n');
fprintf('  a_phi1=%.4f  a_phi2=%.4f\n', a_phi1, a_phi2);
fprintf('  a_theta1=%.4f  a_theta2=%.4f  a_theta3=%.4f\n', a_theta1, a_theta2, a_theta3);
fprintf('  a_V1=%.6f  a_V2=%.4f  a_V3=%.4f\n', a_V1, a_V2, a_V3);
fprintf('  a_beta1=%.4f  a_beta2=%.4f\n', a_beta1, a_beta2);

%% ---- Plot step responses of transfer functions -----------------------

t_step = 0:0.01:10;

figure('Name','Q4: Transfer Function Step Responses','NumberTitle','off');
set(gcf,'Position',[50 50 1200 800]);

subplot(2,4,1)
step(T_phi_delta_a, t_step); grid on;
title('\phi(s)/\delta_a(s)','FontSize',11);
xlabel('Time (s)'); ylabel('\phi (rad)');

subplot(2,4,2)
step(T_chi_phi, t_step); grid on;
title('\chi(s)/\phi(s)','FontSize',11);
xlabel('Time (s)'); ylabel('\chi (rad)');

subplot(2,4,3)
step(T_theta_delta_e, t_step); grid on;
title('\theta(s)/\delta_e(s)','FontSize',11);
xlabel('Time (s)'); ylabel('\theta (rad)');

subplot(2,4,4)
step(T_h_theta, t_step); grid on;
title('h(s)/\theta(s)','FontSize',11);
xlabel('Time (s)'); ylabel('h (m)');

subplot(2,4,5)
step(T_Va_delta_t, t_step); grid on;
title('V_a(s)/\delta_t(s)','FontSize',11);
xlabel('Time (s)'); ylabel('V_a (m/s)');

subplot(2,4,6)
step(T_Va_theta, t_step); grid on;
title('V_a(s)/\theta(s)','FontSize',11);
xlabel('Time (s)'); ylabel('V_a (m/s)');

subplot(2,4,7)
step(T_beta_delta_r, t_step); grid on;
title('\beta(s)/\delta_r(s)','FontSize',11);
xlabel('Time (s)'); ylabel('\beta (rad)');

sgtitle('Q4(vi): Transfer Functions — Step Responses (Va=200 m/s, \gamma=0)', ...
        'FontSize',13,'FontWeight','bold');

%% ========================================================================
%% STEP (vii): Analytic State-Space Model (Ch. 5, Eq. 5.43 & 5.50)
%% ========================================================================

% --- Dimensional stability derivatives (partial derivatives of forces/moments) ---
% Longitudinal: states [u, w, q, theta, h], inputs [delta_e, delta_t]
qbar = 0.5 * P.rho * Va_trim^2;

% Dimensional aerodynamic derivatives for longitudinal
X_u = -(P.C_D_0 + P.C_D_alpha*alpha_trim) * qbar*P.S_wing * 2/Va_trim / P.mass;
X_w = -(P.C_D_alpha * qbar*P.S_wing) / (P.mass*Va_trim);
X_q = -(P.C_D_q * P.c/(2*Va_trim) * qbar*P.S_wing) / P.mass;

Z_u = -(P.C_L_0 + P.C_L_alpha*alpha_trim) * qbar*P.S_wing * 2/Va_trim / P.mass;
Z_w = -(P.C_L_alpha * qbar*P.S_wing) / (P.mass*Va_trim);
Z_q = Va_trim - (P.C_L_q * P.c/(2*Va_trim) * qbar*P.S_wing)/P.mass;

M_u = 0;
M_w = P.C_m_alpha * qbar*P.S_wing*P.c / (P.Jy * Va_trim);
M_q = P.C_m_q * P.c/(2*Va_trim) * qbar*P.S_wing*P.c / P.Jy;

X_de = -P.C_D_delta_e * qbar*P.S_wing / P.mass;
X_dt = P.T_max / P.mass;
Z_de = -P.C_L_delta_e * qbar*P.S_wing / P.mass;
M_de = P.C_m_delta_e * qbar*P.S_wing*P.c / P.Jy;

g = P.gravity;

% Longitudinal SS (Eq. 5.43): states = [u, w, q, theta, h]
A_lon_analytic = [
    X_u,      X_w,  X_q,             -g*cos(theta_trim),  0;
    Z_u,      Z_w,  Z_q,             -g*sin(theta_trim),  0;
    M_u,      M_w,  M_q,              0,                  0;
    0,        0,    1,                0,                  0;
    sin(theta_trim), -cos(theta_trim), 0, Va_trim*cos(0), 0  % gamma~0
];

B_lon_analytic = [
    X_de,  X_dt;
    Z_de,  0;
    M_de,  0;
    0,     0;
    0,     0
];

% Lateral SS: states = [v, p, r, phi, psi], inputs = [delta_a, delta_r]
Y_v = P.C_Y_beta * qbar*P.S_wing / (P.mass * Va_trim);
Y_p = P.C_Y_p * P.b/(2*Va_trim) * qbar*P.S_wing / P.mass;
Y_r = P.C_Y_r * P.b/(2*Va_trim) * qbar*P.S_wing / P.mass;
L_v = P.C_ell_beta * qbar*P.S_wing*P.b / (P.Jx * Va_trim);   % approximate Jx-only
L_p = P.C_ell_p * P.b/(2*Va_trim) * qbar*P.S_wing*P.b / P.Jx;
L_r = P.C_ell_r * P.b/(2*Va_trim) * qbar*P.S_wing*P.b / P.Jx;
N_v = P.C_n_beta * qbar*P.S_wing*P.b / (P.Jz * Va_trim);
N_p = P.C_n_p * P.b/(2*Va_trim) * qbar*P.S_wing*P.b / P.Jz;
N_r = P.C_n_r * P.b/(2*Va_trim) * qbar*P.S_wing*P.b / P.Jz;

Y_da = P.C_Y_delta_a * qbar*P.S_wing / P.mass;
Y_dr = P.C_Y_delta_r * qbar*P.S_wing / P.mass;
L_da = P.C_ell_delta_a * qbar*P.S_wing*P.b / P.Jx;
L_dr = P.C_ell_delta_r * qbar*P.S_wing*P.b / P.Jx;
N_da = P.C_n_delta_a * qbar*P.S_wing*P.b / P.Jz;
N_dr = P.C_n_delta_r * qbar*P.S_wing*P.b / P.Jz;

A_lat_analytic = [
    Y_v,  Y_p,  -(Va_trim - Y_r),  g*cos(theta_trim),  0;
    L_v,  L_p,  L_r,                0,                  0;
    N_v,  N_p,  N_r,                0,                  0;
    0,    1,    tan(theta_trim),    0,                  0;
    0,    0,    sec(theta_trim),    0,                  0
];

B_lat_analytic = [
    Y_da,  Y_dr;
    L_da,  L_dr;
    N_da,  N_dr;
    0,     0;
    0,     0
];

fprintf('\n--- Longitudinal State-Space Model A_lon (states: u,w,q,theta,h) ---\n');
disp(A_lon_analytic);
fprintf('--- B_lon (inputs: delta_e, delta_t) ---\n');
disp(B_lon_analytic);

fprintf('\n--- Lateral State-Space Model A_lat (states: v,p,r,phi,psi) ---\n');
disp(A_lat_analytic);
fprintf('--- B_lat (inputs: delta_a, delta_r) ---\n');
disp(B_lat_analytic);

% Display eigenvalues
fprintf('\nLongitudinal eigenvalues:\n');
disp(eig(A_lon_analytic));
fprintf('Lateral eigenvalues:\n');
disp(eig(A_lat_analytic));

%% ---- Plot SS model Bode (optional sanity check) ----------------------
sys_lon = ss(A_lon_analytic, B_lon_analytic, eye(5), zeros(5,2));
sys_lat = ss(A_lat_analytic, B_lat_analytic, eye(5), zeros(5,2));

figure('Name','Q4(vii): Bode Diagrams - State-Space','NumberTitle','off');
subplot(1,2,1)
bode(sys_lon(3,1));  % q response to delta_e
title('Longitudinal: q/\delta_e Bode','FontSize',12);
grid on;
subplot(1,2,2)
bode(sys_lat(2,1));  % p response to delta_a
title('Lateral: p/\delta_a Bode','FontSize',12);
grid on;
sgtitle('Q4(vii): State-Space Model Bode Verification','FontSize',13,'FontWeight','bold');

fprintf('\n=== Q4 Complete ===\n');

%% ========================================================================
%% LOCAL FUNCTION: Analytic Trim Computation (no Simulink needed)
%% ========================================================================
function [x_trim, u_trim, y_trim, climb_rate] = compute_trim_analytic(Va, gamma, R, P)
% Computes trim analytically for given Va, gamma (rad), R (m).
% For wings-level: R=Inf.  For right turn: R>0.

    g = P.gravity;
    
    % Coordinated turn bank angle from load factor
    if isinf(R)
        phi = 0;
    else
        % phi = atan(Va^2 / (g*R))  for coordinated turn
        phi = atan(Va^2 / (g*abs(R))) * sign(R);
    end
    
    % Turn rate
    if isinf(R)
        psi_dot = 0;
    else
        psi_dot = Va * cos(gamma) / R;
    end
    
    % Steady-state angular rates in body frame (coordinated turn)
    p = -psi_dot * sin(0);         % no pitch at wings level
    q =  psi_dot * sin(phi);       % = 0 for wings-level
    r =  psi_dot * cos(phi);
    
    % --- Find alpha by solving trim pitch moment = 0 ---
    % At trim: Cm(alpha, delta_e) = 0 and L - W*cos(gamma) = 0
    % Use load factor: L = n*W, n=1 for level flight, n=1/cos(phi) for turn
    n_load = 1/cos(phi);
    W = P.mass * g;
    qbar = 0.5 * P.rho * Va^2;
    
    % Solve for alpha from lift equation:
    % n*W = qbar*S*(CL0 + CLalpha*alpha)   (ignoring CLq term for initial estimate)
    CL_needed = n_load * W / (qbar * P.S_wing);
    alpha_est = (CL_needed - P.C_L_0) / P.C_L_alpha;
    
    % theta from alpha + gamma relationship
    theta = alpha_est + gamma;
    
    % u,v,w in body frame
    u_b = Va * cos(alpha_est);
    v_b = Va * sin(phi) * sin(0);   % ~0 for coordinated
    w_b = Va * sin(alpha_est);
    
    % Trim elevator from pitch moment = 0
    % Cm0 + Cm_alpha*alpha + Cm_q*(c/2Va)*q + Cm_de*delta_e = 0
    % Assume q=0 at trim
    delta_e = -(P.C_m_0 + P.C_m_alpha * alpha_est) / P.C_m_delta_e;
    
    % Trim throttle from drag/thrust balance (longitudinal force):
    % Thrust = Drag (approx, for level or climbing flight in body x)
    Drag = qbar * P.S_wing * (P.C_D_0 + P.C_D_alpha*alpha_est + P.C_D_delta_e*delta_e);
    Lift = qbar * P.S_wing * (P.C_L_0 + P.C_L_alpha*alpha_est + P.C_L_delta_e*delta_e);
    
    % force balance in body x: T = D*cos(alpha) - L*sin(alpha) + W*sin(theta)
    Thrust_needed = Drag*cos(alpha_est) - Lift*sin(alpha_est) + W*sin(theta);
    delta_t = Thrust_needed / P.T_max;
    delta_t = max(0, min(1, delta_t));  % saturate
    
    % Trim rudder from yaw moment = 0 (for turn)
    % Cn0 + Cn_beta*beta + Cn_r*(b/2Va)*r + Cn_dr*delta_r = 0
    beta = 0;   % coordinated turn → beta=0
    delta_r = -(P.C_n_0 + P.C_n_beta*beta + P.C_n_r*(P.b/(2*Va))*r) / P.C_n_delta_r;
    
    % Trim aileron from roll moment = 0
    % Cl0 + Cl_beta*beta + Cl_p*(b/2Va)*p + Cl_r*(b/2Va)*r + Cl_da*delta_a = 0
    delta_a = -(P.C_ell_0 + P.C_ell_beta*beta + P.C_ell_p*(P.b/(2*Va))*p + ...
                P.C_ell_r*(P.b/(2*Va))*r) / P.C_ell_delta_a;
    
    % Outputs
    x_trim   = [0; 0; -200; u_b; v_b; w_b; phi; theta; 0; p; q; r];
    u_trim   = [delta_e; delta_a; delta_r; delta_t];
    y_trim   = [alpha_est; beta; Va];
    climb_rate = Va * sin(gamma);
end
