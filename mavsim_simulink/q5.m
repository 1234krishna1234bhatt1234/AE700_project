%% q5.m — Autopilot Design and Simulation (B&M Ch. 6)
clear; clc; close all;

base_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, 'tools'));
addpath(fullfile(base_dir, 'parameters'));

aerosonde_parameters;
simulation_parameters;

%% Trim at wings-level (γ=0), B&M §5.1
Va_trim = 200;    % m/s

[x_trim, u_trim, y_trim, ~] = compute_trim_analytic(Va_trim, 0, Inf, P);
alpha_trim   = y_trim(1);
theta_trim   = x_trim(8);   % state: [pn,pe,pd,u,v,w,phi,theta,psi,p,q,r]
delta_e_trim = u_trim(1);
delta_t_trim = u_trim(4);   % trim throttle ≈ 0.13 for level cruise

g    = P.gravity;
qbar = 0.5*P.rho*Va_trim^2;

%% TF coefficients, B&M Eq. 5.15, 5.22, 5.29–5.33
a_phi1 = -(P.Gamma3*P.C_ell_p + P.Gamma4*P.C_n_p) * (P.rho*Va_trim*P.S_wing*P.b^2/4);
a_phi2 =  (P.Gamma3*P.C_ell_delta_a + P.Gamma4*P.C_n_delta_a) * (P.rho*Va_trim^2*P.S_wing*P.b/2);
a_theta1 = -(P.rho*Va_trim*P.S_wing*P.c^2)/(4*P.Jy)*P.C_m_q;
a_theta2 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy)*P.C_m_alpha;
a_theta3 = (P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy)*P.C_m_delta_e;
C_D_trim = P.C_D_0 + P.C_D_alpha*alpha_trim + P.C_D_delta_e*delta_e_trim;
a_V1 = (P.rho*Va_trim/P.mass)*(P.S_wing*C_D_trim);
a_V2 = P.T_max/P.mass;
a_V3 = g;
a_beta1 = -(P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_beta;
a_beta2 =  (P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_delta_r;

%% Autopilot gains — pole placement, B&M §6.1
zeta = 0.707;

% Roll PD, B&M Eq. 6.3
omega_roll  = 1.2;
kp_roll     = omega_roll^2 / a_phi2;
kd_roll     = (2*zeta*omega_roll - a_phi1) / a_phi2;
delta_a_max = 40*pi/180;

% Course PI — ω_χ = ω_φ/10, B&M Eq. 6.9
omega_course = omega_roll/10;
kp_course    = 2*zeta*omega_course*Va_trim/g;
ki_course    = omega_course^2*Va_trim/g;
phi_max      = 45*pi/180;

% Sideslip PI, B&M Eq. 6.37
omega_beta  = 5;
kp_sideslip = (2*zeta*omega_beta - a_beta1) / a_beta2;
ki_sideslip = omega_beta^2 / a_beta2;
delta_r_max = 40*pi/180;

% Pitch PD, B&M Eq. 6.17–6.18
omega_pitch = 2;
kp_pitch    = (omega_pitch^2 - a_theta2) / a_theta3;
kd_pitch    = (2*zeta*omega_pitch - a_theta1) / a_theta3;
delta_e_max = 40*pi/180;
K_theta_DC  = kp_pitch*a_theta3 / (a_theta2 + kp_pitch*a_theta3);

% Altitude PI — ω_h = ω_θ/10, B&M Eq. 6.23
omega_alt = omega_pitch/10;
kp_alt    = 2*zeta*omega_alt / (K_theta_DC*Va_trim);
ki_alt    = omega_alt^2 / (K_theta_DC*Va_trim);
theta_c_max = 30*pi/180;

% Airspeed via throttle PI, B&M Eq. 6.27–6.28
omega_Vt    = 1;
kp_Va_thr   = (2*zeta*omega_Vt - a_V1) / a_V2;
ki_Va_thr   = omega_Vt^2 / a_V2;

% Airspeed via pitch PI, B&M Eq. 6.31–6.32
omega_Vp    = 0.5;
kp_Va_pit   = (2*zeta*omega_Vp - a_V1) / (-a_V3);
ki_Va_pit   = omega_Vp^2 / (-a_V3);

fprintf('Gains:\n');
fprintf('  Roll:     kp=%.4f  kd=%.4f\n', kp_roll,   kd_roll);
fprintf('  Course:   kp=%.4f  ki=%.4f\n', kp_course, ki_course);
fprintf('  Pitch:    kp=%.4f  kd=%.4f  KthDC=%.4f\n', kp_pitch, kd_pitch, K_theta_DC);
fprintf('  Altitude: kp=%.4f  ki=%.4f\n', kp_alt,    ki_alt);
fprintf('  Va/Thr:   kp=%.4f  ki=%.4f\n', kp_Va_thr, ki_Va_thr);
fprintf('  Va/Pit:   kp=%.4f  ki=%.4f\n', kp_Va_pit, ki_Va_pit);
fprintf('  Trim:     delta_t=%.4f  delta_e=%.4f\n', delta_t_trim, delta_e_trim);

%% Discretise plant TFs (ZOH), B&M §6.1
Ts = SIM.ts_simulation;
Gphi   = c2d(tf(a_phi2,   [1, a_phi1, 0]),          Ts, 'zoh');
Gchi   = c2d(tf(g/Va_trim, [1, 0]),                  Ts, 'zoh');
Gbeta  = c2d(tf(a_beta2,  [1, a_beta1]),             Ts, 'zoh');
Gtheta = c2d(tf(a_theta3, [1, a_theta1, a_theta2]),  Ts, 'zoh');
Gh_th  = c2d(tf(Va_trim,  [1, 0]),                   Ts, 'zoh');
GVa_dt = c2d(tf(a_V2,    [1, a_V1]),                 Ts, 'zoh');
GVa_th = c2d(tf(-a_V3,   [1, a_V1]),                 Ts, 'zoh');

[phi_A,  phi_B,  phi_C,  phi_D]  = ssdata(Gphi);
[chi_A,  chi_B,  chi_C,  chi_D]  = ssdata(Gchi);
[beta_A, beta_B, beta_C, beta_D] = ssdata(Gbeta);
[th_A,   th_B,   th_C,   th_D]   = ssdata(Gtheta);
[h_A,    h_B,    h_C,    h_D]    = ssdata(Gh_th);
[Vt_A,   Vt_B,   Vt_C,   Vt_D]  = ssdata(GVa_dt);
[Vp_A,   Vp_B,   Vp_C,   Vp_D]  = ssdata(GVa_th);

%% Simulation setup
Tend = 120;
t    = 0:Ts:Tend;
N    = length(t);

dist_mag  = 0.2;    % aileron impulse disturbance magnitude (rad), added to delta_a command
dist_time = 30;     % s

chi_cmd_val   = 20*pi/180;
theta_cmd_val = 5*pi/180;
h_cmd_val     = 50;
Va_step       = 10;

%% Lateral — no disturbance
phi_ND     = zeros(1,N);
chi_ND     = zeros(1,N);
phi_cmd_ND = zeros(1,N);
da_ND      = zeros(1,N);
beta_ND    = zeros(1,N);
dr_ND      = zeros(1,N);

phi_x_nd  = zeros(size(phi_A,1),1);
chi_x_nd  = zeros(size(chi_A,1),1);
beta_x_nd = zeros(size(beta_A,1),1);
int_chi_nd = 0;  int_beta_nd = 0;
p_phi_nd   = 0;  phi_prev_nd = 0;

for k = 1:N
    tk = t(k);
    chi_cmd = chi_cmd_val * (tk >= 1.0);

    % Course PI (B&M Eq. 6.7)
    e_chi = chi_cmd - chi_ND(max(1,k-1));
    int_chi_nd = int_chi_nd + Ts*e_chi;
    phi_cmd_ND(k) = sat_fn(kp_course*e_chi + ki_course*int_chi_nd, phi_max, -phi_max);

    % Roll PD (B&M Eq. 6.1)
    e_phi = phi_cmd_ND(k) - phi_ND(max(1,k-1));
    da_ND(k) = sat_fn(kp_roll*e_phi - kd_roll*p_phi_nd, delta_a_max, -delta_a_max);

    phi_x_nd = phi_A*phi_x_nd + phi_B*da_ND(k);
    phi_ND(k) = phi_C*phi_x_nd + phi_D*da_ND(k);
    p_phi_nd  = (phi_ND(k) - phi_prev_nd) / Ts;
    phi_prev_nd = phi_ND(k);

    chi_x_nd = chi_A*chi_x_nd + chi_B*phi_ND(k);
    chi_ND(k) = chi_C*chi_x_nd + chi_D*phi_ND(k);

    % Sideslip PI (B&M Eq. 6.35)
    e_beta = 0 - beta_ND(max(1,k-1));
    int_beta_nd = int_beta_nd + Ts*e_beta;
    dr_ND(k) = sat_fn(kp_sideslip*e_beta + ki_sideslip*int_beta_nd, delta_r_max, -delta_r_max);
    beta_x_nd = beta_A*beta_x_nd + beta_B*dr_ND(k);
    beta_ND(k) = beta_C*beta_x_nd + beta_D*dr_ND(k);
end

%% Lateral — with disturbance
phi_D2     = zeros(1,N);
chi_D2     = zeros(1,N);
phi_cmd_D2 = zeros(1,N);
da_D2      = zeros(1,N);

phi_x_d2  = zeros(size(phi_A,1),1);
chi_x_d2  = zeros(size(chi_A,1),1);
int_chi_d2 = 0;
p_phi_d2   = 0;  phi_prev_d2 = 0;

for k = 1:N
    tk = t(k);
    chi_cmd = chi_cmd_val * (tk >= 1.0);

    e_chi = chi_cmd - chi_D2(max(1,k-1));
    int_chi_d2 = int_chi_d2 + Ts*e_chi;
    phi_cmd_D2(k) = sat_fn(kp_course*e_chi + ki_course*int_chi_d2, phi_max, -phi_max);

    e_phi = phi_cmd_D2(k) - phi_D2(max(1,k-1));
    dist  = dist_mag * (tk >= dist_time && tk < dist_time+Ts);
    da_D2(k) = sat_fn(kp_roll*e_phi - kd_roll*p_phi_d2 + dist, delta_a_max, -delta_a_max);

    phi_x_d2 = phi_A*phi_x_d2 + phi_B*da_D2(k);
    phi_D2(k) = phi_C*phi_x_d2 + phi_D*da_D2(k);
    p_phi_d2  = (phi_D2(k) - phi_prev_d2) / Ts;
    phi_prev_d2 = phi_D2(k);

    chi_x_d2 = chi_A*chi_x_d2 + chi_B*phi_D2(k);
    chi_D2(k) = chi_C*chi_x_d2 + chi_D*phi_D2(k);
end

%% Longitudinal — pitch tracking
theta_PT     = zeros(1,N);
de_PT        = zeros(1,N);
theta_cmd_PT = zeros(1,N);
th_x_PT      = zeros(size(th_A,1),1);
q_theta_PT   = 0;  theta_prev_PT = 0;

for k = 1:N
    tk = t(k);
    theta_cmd_PT(k) = theta_cmd_val * (tk >= 1.0);

    % Pitch PD (B&M Eq. 6.15)
    e_th = theta_cmd_PT(k) - theta_PT(max(1,k-1));
    de_PT(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_theta_PT, delta_e_max, -delta_e_max);

    th_x_PT = th_A*th_x_PT + th_B*de_PT(k);
    theta_PT(k) = th_C*th_x_PT + th_D*de_PT(k);
    q_theta_PT = (theta_PT(k) - theta_prev_PT) / Ts;
    theta_prev_PT = theta_PT(k);
end

%% Longitudinal — altitude hold
% dt_* is throttle PERTURBATION from trim; actual throttle = delta_t_trim + dt_*
theta_AH   = zeros(1,N);
h_AH       = zeros(1,N);
Va_AH      = Va_trim*ones(1,N);
de_AH      = zeros(1,N);
dthrottle_AH      = zeros(1,N);
theta_c_AH = zeros(1,N);

th_x_AH = zeros(size(th_A,1),1);
h_x_AH  = zeros(size(h_A,1),1);
Vt_x_AH = zeros(size(Vt_A,1),1);
int_h_AH = 0;  int_Vt_AH = 0;
q_theta_AH   = 0;  theta_prev_AH = 0;

for k = 1:N
    tk = t(k);
    h_cmd = h_cmd_val * (tk >= 1.0);

    % Altitude PI → θ_c, B&M Eq. 6.20
    e_h = h_cmd - h_AH(max(1,k-1));
    int_h_AH = int_h_AH + Ts*e_h;
    theta_c_AH(k) = sat_fn(kp_alt*e_h + ki_alt*int_h_AH, theta_c_max, -theta_c_max);

    % Pitch PD → δe
    e_th = theta_c_AH(k) - theta_AH(max(1,k-1));
    de_AH(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_theta_AH, delta_e_max, -delta_e_max);

    % Va/throttle PI — dthrottle_AH is perturbation; bounds keep actual δt in [0,1]
    e_Va = Va_trim - Va_AH(max(1,k-1));
    int_Vt_AH = int_Vt_AH + Ts*e_Va;
    dthrottle_AH(k) = sat_fn(kp_Va_thr*e_Va + ki_Va_thr*int_Vt_AH, ...
                      1-delta_t_trim, -delta_t_trim);

    th_x_AH = th_A*th_x_AH + th_B*de_AH(k);
    theta_AH(k) = th_C*th_x_AH + th_D*de_AH(k);
    q_theta_AH = (theta_AH(k) - theta_prev_AH) / Ts;
    theta_prev_AH = theta_AH(k);

    h_x_AH = h_A*h_x_AH + h_B*theta_AH(k);
    h_AH(k) = h_C*h_x_AH + h_D*theta_AH(k);

    Vt_x_AH = Vt_A*Vt_x_AH + Vt_B*dthrottle_AH(k);
    dVa = Vt_C*Vt_x_AH + Vt_D*dthrottle_AH(k);
    Va_AH(k) = Va_trim + dVa;
end

%% Longitudinal — airspeed via throttle
Va_VT    = Va_trim*ones(1,N);
dthrottle_VT    = zeros(1,N);
Vt_x_VT  = zeros(size(Vt_A,1),1);
int_Vt_VT = 0;

for k = 1:N
    tk = t(k);
    Va_cmd = Va_trim + Va_step*(tk >= 1.0);

    e_Va = Va_cmd - Va_VT(max(1,k-1));
    int_Vt_VT = int_Vt_VT + Ts*e_Va;
    dthrottle_VT(k) = sat_fn(kp_Va_thr*e_Va + ki_Va_thr*int_Vt_VT, ...
                      1-delta_t_trim, -delta_t_trim);

    Vt_x_VT = Vt_A*Vt_x_VT + Vt_B*dthrottle_VT(k);
    dVa = Vt_C*Vt_x_VT + Vt_D*dthrottle_VT(k);
    Va_VT(k) = Va_trim + dVa;
end

%% Longitudinal — airspeed via pitch
Va_VP      = zeros(1,N);
theta_VP   = zeros(1,N);
de_VP      = zeros(1,N);
theta_c_VP = zeros(1,N);

th_x_VP = zeros(size(th_A,1),1);
Vp_x_VP = zeros(size(Vp_A,1),1);
int_Vp  = 0;
q_theta_VP  = 0;  theta_prev_VP = 0;

for k = 1:N
    tk = t(k);
    Va_cmd_perturb = Va_step * (tk >= 1.0);

    % Va/pitch PI (B&M Eq. 6.29)
    e_Va = Va_cmd_perturb - Va_VP(max(1,k-1));
    int_Vp = int_Vp + Ts*e_Va;
    theta_c_VP(k) = sat_fn(kp_Va_pit*e_Va + ki_Va_pit*int_Vp, theta_c_max, -theta_c_max);

    e_th = theta_c_VP(k) - theta_VP(max(1,k-1));
    de_VP(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_theta_VP, delta_e_max, -delta_e_max);

    th_x_VP = th_A*th_x_VP + th_B*de_VP(k);
    theta_VP(k) = th_C*th_x_VP + th_D*de_VP(k);
    q_theta_VP = (theta_VP(k) - theta_prev_VP) / Ts;
    theta_prev_VP = theta_VP(k);

    Vp_x_VP = Vp_A*Vp_x_VP + Vp_B*theta_VP(k);
    Va_VP(k) = Vp_C*Vp_x_VP + Vp_D*theta_VP(k);
end

%% Precompute actual throttle (trim + perturbation) for plots
throttle_actual_AH = delta_t_trim + dthrottle_AH;
throttle_actual_VT = delta_t_trim + dthrottle_VT;

%% Plots — Lateral

figure('Name','Q5-LAT: Roll — No Disturbance','NumberTitle','off');
set(gcf,'Position',[50 600 900 400]);
subplot(1,2,1)
plot(t, phi_ND*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, phi_cmd_ND*180/pi, 'r--', 'LineWidth',1.5);
legend('\phi actual','\phi commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\phi (deg)','FontSize',11);
title('Roll Angle — No Disturbance','FontSize',12); grid on;
subplot(1,2,2)
plot(t, da_ND*180/pi, 'b-', 'LineWidth',1.5);
yline( delta_a_max*180/pi,'r--','Max','FontSize',9);
yline(-delta_a_max*180/pi,'r--','Min','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_a (deg)','FontSize',11);
title('Aileron — No Disturbance','FontSize',12); grid on;
sgtitle('Q5 Lateral: Roll Angle Tracking (No Disturbance)','FontSize',13,'FontWeight','bold');

figure('Name','Q5-LAT: Roll — With Disturbance','NumberTitle','off');
set(gcf,'Position',[50 100 900 400]);
subplot(1,2,1)
plot(t, phi_D2*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, phi_cmd_D2*180/pi, 'r--', 'LineWidth',1.5);
xline(dist_time,'k:','Disturbance','LabelVerticalAlignment','bottom','FontSize',9);
legend('\phi actual','\phi commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\phi (deg)','FontSize',11);
title('Roll Angle — With Disturbance (step = 0.2 rad)','FontSize',12); grid on;
subplot(1,2,2)
plot(t, da_D2*180/pi, 'b-', 'LineWidth',1.5);
xline(dist_time,'k:','FontSize',9);
yline( delta_a_max*180/pi,'r--','Max','FontSize',9);
yline(-delta_a_max*180/pi,'r--','Min','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_a (deg)','FontSize',11);
title('Aileron — With Disturbance','FontSize',12); grid on;
sgtitle('Q5 Lateral: Roll Angle Tracking (With Disturbance)','FontSize',13,'FontWeight','bold');

figure('Name','Q5-LAT: Course Angle','NumberTitle','off');
set(gcf,'Position',[1000 600 700 380]);
plot(t, chi_ND*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, chi_D2*180/pi, 'm-', 'LineWidth',1.5);
yline(chi_cmd_val*180/pi, 'r--', 'Commanded','FontSize',10);
legend('\chi (no dist)','\chi (with dist)','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\chi (deg)','FontSize',11);
title('Q5 Lateral: Course Angle vs. Commanded','FontSize',13,'FontWeight','bold'); grid on;

figure('Name','Q5-LAT: Sideslip Angle','NumberTitle','off');
set(gcf,'Position',[1000 150 700 380]);
subplot(1,2,1)
plot(t, beta_ND*180/pi, 'b-', 'LineWidth',1.5); hold on;
yline(0,'r--','Commanded = 0','FontSize',9);
legend('\beta actual','\beta commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\beta (deg)','FontSize',11);
title('Sideslip Angle','FontSize',12); grid on;
subplot(1,2,2)
plot(t, dr_ND*180/pi, 'b-', 'LineWidth',1.5);
yline( delta_r_max*180/pi,'r--','FontSize',9);
yline(-delta_r_max*180/pi,'r--','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_r (deg)','FontSize',11);
title('Rudder Deflection','FontSize',12); grid on;
sgtitle('Q5 Lateral: Sideslip & Rudder','FontSize',13,'FontWeight','bold');

%% Plots — Longitudinal

figure('Name','Q5-LON: Pitch Tracking','NumberTitle','off');
set(gcf,'Position',[50 600 900 380]);
subplot(1,2,1)
plot(t, theta_PT*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, theta_cmd_PT*180/pi, 'r--', 'LineWidth',1.5);
legend('\theta actual','\theta commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\theta (deg)','FontSize',11);
title('Pitch Angle Tracking','FontSize',12); grid on;
subplot(1,2,2)
plot(t, de_PT*180/pi, 'b-', 'LineWidth',1.5);
yline( delta_e_max*180/pi,'r--','FontSize',9);
yline(-delta_e_max*180/pi,'r--','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_e (deg)','FontSize',11);
title('Elevator Deflection (Pitch Loop)','FontSize',12); grid on;
sgtitle('Q5 Longitudinal: Pitch Angle vs. Commanded','FontSize',13,'FontWeight','bold');

figure('Name','Q5-LON: Altitude Hold','NumberTitle','off');
set(gcf,'Position',[50 100 1100 380]);
subplot(1,3,1)
plot(t, h_AH, 'b-', 'LineWidth',1.5); hold on;
plot(t, h_cmd_val*(t>=1), 'r--', 'LineWidth',1.5);
legend('h actual','h commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('h (m)','FontSize',11);
title('Altitude Hold','FontSize',12); grid on;
subplot(1,3,2)
plot(t, theta_c_AH*180/pi, 'b-', 'LineWidth',1.5);
xlabel('Time (s)','FontSize',11); ylabel('\theta_c (deg)','FontSize',11);
title('Pitch Command (Altitude Loop)','FontSize',12); grid on;
subplot(1,3,3)
plot(t, de_AH*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, throttle_actual_AH, 'g-', 'LineWidth',1.5);
yline( delta_e_max*180/pi,'r:','FontSize',9);
yline(-delta_e_max*180/pi,'r:','FontSize',9);
yline(1,'k:','FontSize',9); yline(0,'k:','FontSize',9);
legend('\delta_e (deg)','\delta_t actual','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('Control','FontSize',11);
title('Elevator & Throttle (actual)','FontSize',12); grid on;
sgtitle('Q5 Longitudinal: Altitude Hold','FontSize',13,'FontWeight','bold');

figure('Name','Q5-LON: Airspeed via Throttle','NumberTitle','off');
set(gcf,'Position',[1050 600 900 380]);
subplot(1,2,1)
plot(t, Va_VT, 'b-', 'LineWidth',1.5); hold on;
plot(t, Va_trim + Va_step*(t>=1), 'r--', 'LineWidth',1.5);
legend('V_a actual','V_a commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('V_a (m/s)','FontSize',11);
title('Airspeed via Throttle','FontSize',12); grid on;
subplot(1,2,2)
plot(t, throttle_actual_VT, 'b-', 'LineWidth',1.5);
yline(1,'r:','Max','FontSize',9); yline(0,'r:','Min','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_t (actual)','FontSize',11);
title('Throttle Setting (actual)','FontSize',12); grid on;
sgtitle('Q5 Longitudinal: Airspeed Control via Throttle','FontSize',13,'FontWeight','bold');

figure('Name','Q5-LON: Airspeed via Pitch','NumberTitle','off');
set(gcf,'Position',[1050 100 900 380]);
subplot(1,2,1)
plot(t, Va_VP, 'b-', 'LineWidth',1.5); hold on;
plot(t, Va_step*(t>=1), 'r--', 'LineWidth',1.5);
legend('\DeltaV_a actual','\DeltaV_a commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\DeltaV_a (m/s)','FontSize',11);
title('Airspeed via Pitch (perturbation)','FontSize',12); grid on;
subplot(1,2,2)
plot(t, de_VP*180/pi, 'b-', 'LineWidth',1.5);
yline( delta_e_max*180/pi,'r--','FontSize',9);
yline(-delta_e_max*180/pi,'r--','FontSize',9);
xlabel('Time (s)','FontSize',11); ylabel('\delta_e (deg)','FontSize',11);
title('Elevator (Airspeed via Pitch)','FontSize',12); grid on;
sgtitle('Q5 Longitudinal: Airspeed Control via Pitch Command','FontSize',13,'FontWeight','bold');

fprintf('\n=== Q5 Autopilot Simulation Complete. Check Figures 1–8. ===\n');

%% Local functions

function out = sat_fn(in, hi, lo)
    out = max(lo, min(hi, in));
end

function [x_trim, u_trim, y_trim, climb_rate] = compute_trim_analytic(Va, gamma, R, P)
    % Analytic trim for given Va (m/s), gamma (rad), R (m). R=Inf → wings-level.
    % B&M §5.1, Eqs. 5.4–5.8.
    g = P.gravity;
    if isinf(R)
        phi = 0;  psi_dot = 0;
    else
        phi     = atan(Va^2/(g*abs(R))) * sign(R);
        psi_dot = Va*cos(gamma)/R;
    end
    p = 0;
    q = psi_dot*sin(phi);
    r = psi_dot*cos(phi);

    n_load    = 1/cos(phi);
    W         = P.mass*g;
    qbar      = 0.5*P.rho*Va^2;
    CL_needed = n_load*W/(qbar*P.S_wing);
    alpha_est = (CL_needed - P.C_L_0)/P.C_L_alpha;
    theta     = alpha_est + gamma;

    u_b = Va*cos(alpha_est);
    v_b = 0;
    w_b = Va*sin(alpha_est);

    % Elevator from pitch moment trim (B&M Eq. 5.6)
    delta_e = -(P.C_m_0 + P.C_m_alpha*alpha_est)/P.C_m_delta_e;

    % Throttle from longitudinal force trim (B&M Eq. 5.4)
    Drag = qbar*P.S_wing*(P.C_D_0 + P.C_D_alpha*alpha_est + P.C_D_delta_e*delta_e);
    Lift = qbar*P.S_wing*(P.C_L_0 + P.C_L_alpha*alpha_est + P.C_L_delta_e*delta_e);
    T_needed = Drag*cos(alpha_est) - Lift*sin(alpha_est) + W*sin(theta);
    delta_t  = max(0, min(1, T_needed/P.T_max));

    % Rudder from yaw moment trim (B&M Eq. 5.8)
    beta    = 0;
    delta_r = -(P.C_n_0 + P.C_n_beta*beta + P.C_n_r*(P.b/(2*Va))*r) / P.C_n_delta_r;

    % Aileron from roll moment trim (B&M Eq. 5.7)
    delta_a = -(P.C_ell_0 + P.C_ell_beta*beta ...
                + P.C_ell_p*(P.b/(2*Va))*p ...
                + P.C_ell_r*(P.b/(2*Va))*r) / P.C_ell_delta_a;

    x_trim    = [0; 0; -200; u_b; v_b; w_b; phi; theta; 0; p; q; r];
    u_trim    = [delta_e; delta_a; delta_r; delta_t];
    y_trim    = [alpha_est; beta; Va];
    climb_rate = Va*sin(gamma);
end
