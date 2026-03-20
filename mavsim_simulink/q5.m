%% ========================================================================
%% q5.m  —  Autopilot Design & Simulation (Q5 of AE700 project)
%% ========================================================================
%% Simulates lateral and longitudinal autopilots.
%% Produces all required plots (with/without disturbance).
%% ========================================================================
clear; clc; close all;

%% ---- Resolve paths robustly -------------------------------------------
base_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(base_dir, 'tools'));       % add tools BEFORE aerosonde_parameters
addpath(fullfile(base_dir, 'parameters'));

aerosonde_parameters;
simulation_parameters;

%% ---- Compute trim (wings-level, gamma=0) -------------------------------
Va_trim = 200;    % m/s (from aerosonde_parameters)

[x_trim, u_trim, y_trim, ~] = compute_trim_analytic(Va_trim, 0, Inf, P);
alpha_trim   = y_trim(1);
theta_trim   = x_trim(8);
delta_e_trim = u_trim(1);

g    = P.gravity;
qbar = 0.5*P.rho*Va_trim^2;

%% ---- Transfer function coefficients ------------------------------------
a_phi1 = -(P.Gamma3*P.C_ell_p + P.Gamma4*P.C_n_p) ...
          * (P.rho*Va_trim*P.S_wing*P.b^2/4);
a_phi2 =  (P.Gamma3*P.C_ell_delta_a + P.Gamma4*P.C_n_delta_a) ...
          * (P.rho*Va_trim^2*P.S_wing*P.b/2);
a_theta1 = -(P.rho*Va_trim*P.S_wing*P.c^2)/(4*P.Jy)*P.C_m_q;
a_theta2 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy)*P.C_m_alpha;
a_theta3 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy)*P.C_m_delta_e;
C_D_trim = P.C_D_0 + P.C_D_alpha*alpha_trim + P.C_D_delta_e*delta_e_trim;
a_V1 = (P.rho*Va_trim/P.mass)*(P.S_wing*C_D_trim);
a_V2 = P.T_max/P.mass;
a_V3 = g;
a_beta1 = -(P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_beta;
a_beta2 =  (P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_delta_r;

%% ---- Autopilot gains (pole placement) ----------------------------------
zeta = 0.707;

% Roll (PD)
omega_roll  = 10;
kp_roll     = omega_roll^2 / a_phi2;
kd_roll     = (2*zeta*omega_roll - a_phi1) / a_phi2;
delta_a_max = 40*pi/180;

% Course (PI), ~10x slower
omega_course = omega_roll/10;
kp_course    = 2*zeta*omega_course*Va_trim/g;
ki_course    = omega_course^2*Va_trim/g;
phi_max      = 45*pi/180;

% Sideslip (PI)
omega_beta  = 5;
kp_sideslip = (2*zeta*omega_beta - a_beta1) / a_beta2;
ki_sideslip = omega_beta^2 / a_beta2;
delta_r_max = 40*pi/180;

% Yaw damper (P)
kp_yaw_damper = 0.5;

% Pitch (PD)
omega_pitch = 10;
kp_pitch    = (omega_pitch^2 - a_theta2) / a_theta3;
kd_pitch    = (2*zeta*omega_pitch - a_theta1) / a_theta3;
delta_e_max = 40*pi/180;
K_theta_DC  = kp_pitch*a_theta3 / (a_theta2 + kp_pitch*a_theta3);

% Altitude (PI)
omega_alt = omega_pitch/10;
kp_alt    = 2*zeta*omega_alt / (K_theta_DC*Va_trim);
ki_alt    = omega_alt^2 / (K_theta_DC*Va_trim);
theta_c_max = 30*pi/180;

% Airspeed via throttle (PI)
omega_Vt    = 2;
kp_Va_thr   = (2*zeta*omega_Vt - a_V1) / a_V2;
ki_Va_thr   = omega_Vt^2 / a_V2;

% Airspeed via pitch (PI)
omega_Vp    = 1;
kp_Va_pit   = (2*zeta*omega_Vp - a_V1) / (-a_V3);
ki_Va_pit   = omega_Vp^2 / (-a_V3);

fprintf('Gains:\n');
fprintf('  Roll:     kp=%.4f  kd=%.4f\n', kp_roll,   kd_roll);
fprintf('  Course:   kp=%.4f  ki=%.4f\n', kp_course, ki_course);
fprintf('  Pitch:    kp=%.4f  kd=%.4f  KthDC=%.4f\n', kp_pitch, kd_pitch, K_theta_DC);
fprintf('  Altitude: kp=%.4f  ki=%.4f\n', kp_alt,    ki_alt);
fprintf('  Va/Thr:   kp=%.4f  ki=%.4f\n', kp_Va_thr, ki_Va_thr);
fprintf('  Va/Pit:   kp=%.4f  ki=%.4f\n', kp_Va_pit, ki_Va_pit);

%% ---- Discretise plant transfer functions --------------------------------
Ts = SIM.ts_simulation;   % 0.02 s
Gphi  = c2d(tf(a_phi2,   [1, a_phi1, 0]),           Ts, 'zoh');
Gchi  = c2d(tf(g/Va_trim, [1, 0]),                   Ts, 'zoh');
Gbeta = c2d(tf(a_beta2,  [1, a_beta1]),              Ts, 'zoh');
Gtheta = c2d(tf(a_theta3, [1, a_theta1, a_theta2]),  Ts, 'zoh');
Gh_th  = c2d(tf(Va_trim,  [1, 0]),                   Ts, 'zoh');
GVa_dt = c2d(tf(a_V2,    [1, a_V1]),                 Ts, 'zoh');
GVa_th = c2d(tf(-a_V3,   [1, a_V1]),                 Ts, 'zoh');

[phi_A, phi_B, phi_C, phi_D]   = ssdata(Gphi);
[chi_A, chi_B, chi_C, chi_D]   = ssdata(Gchi);
[beta_A, beta_B, beta_C, beta_D] = ssdata(Gbeta);
[th_A, th_B, th_C, th_D]       = ssdata(Gtheta);
[h_A,  h_B,  h_C,  h_D]        = ssdata(Gh_th);
[Vt_A, Vt_B, Vt_C, Vt_D]       = ssdata(GVa_dt);
[Vp_A, Vp_B, Vp_C, Vp_D]       = ssdata(GVa_th);

%% ========================================================================
%%  SIMULATION PARAMETERS
%% ========================================================================
Tend = 60;
t    = 0:Ts:Tend;
N    = length(t);

dist_mag  = 0.2;    % rad step disturbance
dist_time = 15;     % seconds

chi_cmd_val   = 20*pi/180;  % 20 deg course step
theta_cmd_val = 5*pi/180;   % 5 deg pitch step
h_cmd_val     = 50;          % 50 m altitude step
Va_step       = 10;          % m/s airspeed step above trim

%% ========================================================================
%%  LATERAL — NO DISTURBANCE
%% ========================================================================
phi_ND    = zeros(1,N);
chi_ND    = zeros(1,N);
phi_cmd_ND = zeros(1,N);
da_ND     = zeros(1,N);
beta_ND   = zeros(1,N);
dr_ND     = zeros(1,N);

phi_x_nd  = zeros(size(phi_A,1),1);
chi_x_nd  = zeros(size(chi_A,1),1);
beta_x_nd = zeros(size(beta_A,1),1);
int_chi_nd  = 0;
int_beta_nd = 0;

for k = 1:N
    tk = t(k);
    chi_cmd = chi_cmd_val * (tk >= 1.0);

    % --- Course PI → phi command ---
    e_chi = chi_cmd - chi_ND(max(1,k-1));
    int_chi_nd = int_chi_nd + Ts*e_chi;
    phi_cmd_ND(k) = sat_fn(kp_course*e_chi + ki_course*int_chi_nd, phi_max, -phi_max);

    % --- Roll PD → delta_a ---
    phi_prev = phi_ND(max(1,k-1));
    % FIX: use max(1,k-1) and max(1,k-2) so index never reaches 0
    p_nd = (phi_ND(max(1,k-1)) - phi_ND(max(1,k-2))) / Ts * (k>1);
    e_phi = phi_cmd_ND(k) - phi_prev;
    da_ND(k) = sat_fn(kp_roll*e_phi - kd_roll*p_nd, delta_a_max, -delta_a_max);

    % --- Roll plant ---
    phi_x_nd = phi_A*phi_x_nd + phi_B*da_ND(k);
    phi_ND(k) = phi_C*phi_x_nd + phi_D*da_ND(k);

    % --- Course plant ---
    chi_x_nd = chi_A*chi_x_nd + chi_B*phi_ND(k);
    chi_ND(k) = chi_C*chi_x_nd + chi_D*phi_ND(k);

    % --- Sideslip PI → delta_r ---
    beta_prev = beta_ND(max(1,k-1));
    e_beta = 0 - beta_prev;
    int_beta_nd = int_beta_nd + Ts*e_beta;
    dr_ND(k) = sat_fn(kp_sideslip*e_beta + ki_sideslip*int_beta_nd, delta_r_max, -delta_r_max);
    beta_x_nd = beta_A*beta_x_nd + beta_B*dr_ND(k);
    beta_ND(k) = beta_C*beta_x_nd + beta_D*dr_ND(k);
end

%% ========================================================================
%%  LATERAL — WITH DISTURBANCE
%% ========================================================================
phi_D2    = zeros(1,N);
chi_D2    = zeros(1,N);
phi_cmd_D2 = zeros(1,N);
da_D2     = zeros(1,N);

phi_x_d2  = zeros(size(phi_A,1),1);
chi_x_d2  = zeros(size(chi_A,1),1);
int_chi_d2 = 0;

for k = 1:N
    tk = t(k);
    chi_cmd = chi_cmd_val * (tk >= 1.0);

    e_chi = chi_cmd - chi_D2(max(1,k-1));
    int_chi_d2 = int_chi_d2 + Ts*e_chi;
    phi_cmd_D2(k) = sat_fn(kp_course*e_chi + ki_course*int_chi_d2, phi_max, -phi_max);

    phi_prev = phi_D2(max(1,k-1));
    % FIX: safe index for rate
    p_d2 = (phi_D2(max(1,k-1)) - phi_D2(max(1,k-2))) / Ts * (k>1);
    e_phi = phi_cmd_D2(k) - phi_prev;

    % Step disturbance
    dist = dist_mag * (tk >= dist_time && tk < dist_time+Ts);

    da_D2(k) = sat_fn(kp_roll*e_phi - kd_roll*p_d2 + dist, delta_a_max, -delta_a_max);

    phi_x_d2 = phi_A*phi_x_d2 + phi_B*da_D2(k);
    phi_D2(k) = phi_C*phi_x_d2 + phi_D*da_D2(k);

    chi_x_d2 = chi_A*chi_x_d2 + chi_B*phi_D2(k);
    chi_D2(k) = chi_C*chi_x_d2 + chi_D*phi_D2(k);
end

%% ========================================================================
%%  LONGITUDINAL — PITCH TRACKING
%% ========================================================================
theta_PT   = zeros(1,N);
de_PT      = zeros(1,N);
theta_cmd_PT = zeros(1,N);
th_x_PT    = zeros(size(th_A,1),1);

for k = 1:N
    tk = t(k);
    theta_cmd_PT(k) = theta_cmd_val * (tk >= 1.0);

    theta_prev = theta_PT(max(1,k-1));
    % FIX: safe index for pitch rate
    q_PT = (theta_PT(max(1,k-1)) - theta_PT(max(1,k-2))) / Ts * (k>1);
    e_th = theta_cmd_PT(k) - theta_prev;
    de_PT(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_PT, delta_e_max, -delta_e_max);

    th_x_PT = th_A*th_x_PT + th_B*de_PT(k);
    theta_PT(k) = th_C*th_x_PT + th_D*de_PT(k);
end

%% ========================================================================
%%  LONGITUDINAL — ALTITUDE HOLD (throttle → Va, pitch → h)
%% ========================================================================
theta_AH  = zeros(1,N);
h_AH      = zeros(1,N);
Va_AH     = Va_trim*ones(1,N);
de_AH     = zeros(1,N);
dt_AH     = zeros(1,N);
theta_c_AH = zeros(1,N);

th_x_AH = zeros(size(th_A,1),1);
h_x_AH  = zeros(size(h_A,1),1);
Vt_x_AH = zeros(size(Vt_A,1),1);
int_h_AH  = 0;
int_Vt_AH = 0;

for k = 1:N
    tk = t(k);
    h_cmd = h_cmd_val * (tk >= 1.0);

    % Altitude PI → theta_c
    e_h = h_cmd - h_AH(max(1,k-1));
    int_h_AH = int_h_AH + Ts*e_h;
    theta_c_AH(k) = sat_fn(kp_alt*e_h + ki_alt*int_h_AH, theta_c_max, -theta_c_max);

    % Pitch PD → delta_e
    theta_prev = theta_AH(max(1,k-1));
    % FIX: safe index for pitch rate
    q_AH = (theta_AH(max(1,k-1)) - theta_AH(max(1,k-2))) / Ts * (k>1);
    e_th = theta_c_AH(k) - theta_prev;
    de_AH(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_AH, delta_e_max, -delta_e_max);

    % Airspeed/throttle PI
    e_Va = Va_trim - Va_AH(max(1,k-1));
    int_Vt_AH = int_Vt_AH + Ts*e_Va;
    dt_AH(k) = sat_fn(kp_Va_thr*e_Va + ki_Va_thr*int_Vt_AH, 1, 0);

    % Update plants
    th_x_AH = th_A*th_x_AH + th_B*de_AH(k);
    theta_AH(k) = th_C*th_x_AH + th_D*de_AH(k);

    h_x_AH = h_A*h_x_AH + h_B*theta_AH(k);
    h_AH(k) = h_C*h_x_AH + h_D*theta_AH(k);

    Vt_x_AH = Vt_A*Vt_x_AH + Vt_B*dt_AH(k);
    dVa = Vt_C*Vt_x_AH + Vt_D*dt_AH(k);
    Va_AH(k) = Va_trim + dVa;
end

%% ========================================================================
%%  LONGITUDINAL — AIRSPEED VIA THROTTLE (constant altitude, Va step)
%% ========================================================================
Va_VT  = Va_trim*ones(1,N);
dt_VT  = zeros(1,N);
Vt_x_VT = zeros(size(Vt_A,1),1);
int_Vt_VT = 0;

for k = 1:N
    tk = t(k);
    Va_cmd = Va_trim + Va_step*(tk >= 1.0);

    e_Va = Va_cmd - Va_VT(max(1,k-1));
    int_Vt_VT = int_Vt_VT + Ts*e_Va;
    dt_VT(k) = sat_fn(kp_Va_thr*e_Va + ki_Va_thr*int_Vt_VT, 1, 0);

    Vt_x_VT = Vt_A*Vt_x_VT + Vt_B*dt_VT(k);
    dVa = Vt_C*Vt_x_VT + Vt_D*dt_VT(k);
    Va_VT(k) = Va_trim + dVa;
end

%% ========================================================================
%%  LONGITUDINAL — AIRSPEED VIA PITCH (Va step using theta command)
%% ========================================================================
Va_VP     = zeros(1,N);     % perturbation around trim
theta_VP  = zeros(1,N);
de_VP     = zeros(1,N);
theta_c_VP = zeros(1,N);

th_x_VP  = zeros(size(th_A,1),1);
Vp_x_VP  = zeros(size(Vp_A,1),1);
int_Vp   = 0;

for k = 1:N
    tk = t(k);
    Va_cmd_perturb = Va_step * (tk >= 1.0);   % perturbation from trim

    e_Va = Va_cmd_perturb - Va_VP(max(1,k-1));
    int_Vp = int_Vp + Ts*e_Va;
    theta_c_VP(k) = sat_fn(kp_Va_pit*e_Va + ki_Va_pit*int_Vp, theta_c_max, -theta_c_max);

    theta_prev = theta_VP(max(1,k-1));
    % FIX: safe index for pitch rate
    q_VP = (theta_VP(max(1,k-1)) - theta_VP(max(1,k-2))) / Ts * (k>1);
    e_th = theta_c_VP(k) - theta_prev;
    de_VP(k) = sat_fn(kp_pitch*e_th - kd_pitch*q_VP, delta_e_max, -delta_e_max);

    th_x_VP = th_A*th_x_VP + th_B*de_VP(k);
    theta_VP(k) = th_C*th_x_VP + th_D*de_VP(k);

    Vp_x_VP = Vp_A*Vp_x_VP + Vp_B*theta_VP(k);
    Va_VP(k) = Vp_C*Vp_x_VP + Vp_D*theta_VP(k);
end

%% ========================================================================
%%  PLOTS — LATERAL
%% ========================================================================

%--- Figure 1: Roll — no disturbance ---
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

%--- Figure 2: Roll — with disturbance ---
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

%--- Figure 3: Course angle ---
figure('Name','Q5-LAT: Course Angle','NumberTitle','off');
set(gcf,'Position',[1000 600 700 380]);

plot(t, chi_ND*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, chi_D2*180/pi, 'm-', 'LineWidth',1.5);
yline(chi_cmd_val*180/pi, 'r--', 'Commanded','FontSize',10);
legend('\chi (no dist)','\chi (with dist)','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\chi (deg)','FontSize',11);
title('Q5 Lateral: Course Angle vs. Commanded','FontSize',13,'FontWeight','bold'); grid on;

%--- Figure 4: Sideslip angle ---
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

%% ========================================================================
%%  PLOTS — LONGITUDINAL
%% ========================================================================

%--- Figure 5: Pitch tracking ---
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

%--- Figure 6: Altitude hold ---
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
plot(t, dt_AH, 'g-', 'LineWidth',1.5);
yline( delta_e_max*180/pi,'r:','FontSize',9);
yline(-delta_e_max*180/pi,'r:','FontSize',9);
legend('\delta_e (deg)','\delta_t (throttle)','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('Control','FontSize',11);
title('Elevator & Throttle','FontSize',12); grid on;

sgtitle('Q5 Longitudinal: Altitude Hold','FontSize',13,'FontWeight','bold');

%--- Figure 7: Airspeed via throttle ---
figure('Name','Q5-LON: Airspeed via Throttle','NumberTitle','off');
set(gcf,'Position',[1050 600 900 380]);

subplot(1,2,1)
plot(t, Va_VT, 'b-', 'LineWidth',1.5); hold on;
plot(t, Va_trim + Va_step*(t>=1), 'r--', 'LineWidth',1.5);
legend('V_a actual','V_a commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('V_a (m/s)','FontSize',11);
title('Airspeed via Throttle','FontSize',12); grid on;

subplot(1,2,2)
plot(t, dt_VT, 'b-', 'LineWidth',1.5);
ylim([0 1.05]);
xlabel('Time (s)','FontSize',11); ylabel('\delta_t','FontSize',11);
title('Throttle Setting','FontSize',12); grid on;

sgtitle('Q5 Longitudinal: Airspeed Control via Throttle','FontSize',13,'FontWeight','bold');

%--- Figure 8: Airspeed via pitch ---
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

%% ========================================================================
%%  LOCAL FUNCTIONS
%% ========================================================================

function out = sat_fn(in, hi, lo)
    out = max(lo, min(hi, in));
end

function [x_trim, u_trim, y_trim, climb_rate] = compute_trim_analytic(Va, gamma, R, P)
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

    delta_e = -(P.C_m_0 + P.C_m_alpha*alpha_est)/P.C_m_delta_e;
    Drag = qbar*P.S_wing*(P.C_D_0 + P.C_D_alpha*alpha_est + P.C_D_delta_e*delta_e);
    Lift = qbar*P.S_wing*(P.C_L_0 + P.C_L_alpha*alpha_est + P.C_L_delta_e*delta_e);
    T_needed = Drag*cos(alpha_est) - Lift*sin(alpha_est) + W*sin(theta);
    delta_t  = max(0, min(1, T_needed/P.T_max));

    beta    = 0;
    delta_r = -(P.C_n_0 + P.C_n_beta*beta + P.C_n_r*(P.b/(2*Va))*r) / P.C_n_delta_r;
    delta_a = -(P.C_ell_0 + P.C_ell_beta*beta ...
                + P.C_ell_p*(P.b/(2*Va))*p ...
                + P.C_ell_r*(P.b/(2*Va))*r) / P.C_ell_delta_a;

    x_trim    = [0; 0; -200; u_b; v_b; w_b; phi; theta; 0; p; q; r];
    u_trim    = [delta_e; delta_a; delta_r; delta_t];
    y_trim    = [alpha_est; beta; Va];
    climb_rate = Va*sin(gamma);
end
