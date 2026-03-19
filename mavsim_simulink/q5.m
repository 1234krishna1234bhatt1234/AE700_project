%% ========================================================================
%% Q5_autopilot_simulation.m
%% ========================================================================
%% Simulates lateral and longitudinal autopilots using linear state-space
%% models derived in Q4.  Produces all required plots from the problem.
%% ========================================================================
clear; clc; close all;

%% ---- Load parameters & run Q4 to get trim + TF coefficients ----------
addpath('parameters');
addpath('tools');
aerosonde_parameters;
simulation_parameters;

Va_trim   = 200;          % m/s
gamma_0   = 0;
R_inf     = Inf;

% Compute trim analytically
[x_trim, u_trim, y_trim, ~] = compute_trim_analytic(Va_trim, gamma_0, R_inf, P);
alpha_trim   = y_trim(1);
theta_trim   = x_trim(8);
delta_e_trim = u_trim(1);
delta_t_trim = u_trim(4);

qbar = 0.5*P.rho*Va_trim^2;
g    = P.gravity;

%% ---- Compute TF coefficients -----------------------------------------
a_phi1 = -(P.Gamma3*P.C_ell_p + P.Gamma4*P.C_n_p) ...
          * (P.rho*Va_trim*P.S_wing*P.b^2/4);
a_phi2 =  (P.Gamma3*P.C_ell_delta_a + P.Gamma4*P.C_n_delta_a) ...
          * (P.rho*Va_trim^2*P.S_wing*P.b/2);
a_theta1 = -(P.rho*Va_trim*P.S_wing*P.c^2)/(4*P.Jy) * P.C_m_q;
a_theta2 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy) * P.C_m_alpha;
a_theta3 = -(P.rho*Va_trim^2*P.S_wing*P.c)/(2*P.Jy) * P.C_m_delta_e;
C_D_trim = P.C_D_0 + P.C_D_alpha*alpha_trim + P.C_D_delta_e*delta_e_trim;
a_V1 = (P.rho*Va_trim/P.mass)*(P.S_wing*C_D_trim);
a_V2 = P.T_max/P.mass;
a_V3 = g;
a_beta1 = -(P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_beta;
a_beta2 =  (P.rho*Va_trim*P.S_wing)/(2*P.mass)*P.C_Y_delta_r;

%% ---- Design Autopilot Gains (pole placement) -------------------------
zeta = 0.707;

% Roll loop (PD)
omega_roll   = 10;
kp_roll      = omega_roll^2 / a_phi2;
kd_roll      = (2*zeta*omega_roll - a_phi1) / a_phi2;
delta_a_max  = 40*pi/180;

% Course loop (PI), ~10x slower
omega_course = omega_roll/10;
kp_course    = 2*zeta*omega_course*Va_trim/g;
ki_course    = omega_course^2*Va_trim/g;
phi_max      = 45*pi/180;

% Yaw damper (P)
kp_yaw_damper = 0.5;
delta_r_max   = 40*pi/180;

% Sideslip PI
omega_beta   = 5;
kp_sideslip  = (2*zeta*omega_beta - a_beta1) / a_beta2;
ki_sideslip  = omega_beta^2 / a_beta2;

% Pitch loop (PD)
omega_pitch  = 10;
kp_pitch     = (omega_pitch^2 - a_theta2) / a_theta3;
kd_pitch     = (2*zeta*omega_pitch - a_theta1) / a_theta3;
delta_e_max  = 40*pi/180;
K_theta_DC   = kp_pitch*a_theta3 / (a_theta2 + kp_pitch*a_theta3);

% Altitude loop (PI), ~10x slower
omega_alt    = omega_pitch/10;
kp_alt       = 2*zeta*omega_alt / (K_theta_DC*Va_trim);
ki_alt       = omega_alt^2 / (K_theta_DC*Va_trim);
theta_c_max  = 30*pi/180;

% Va/throttle (PI)
omega_Vt     = 2;
kp_Va_thr    = (2*zeta*omega_Vt - a_V1) / a_V2;
ki_Va_thr    = omega_Vt^2 / a_V2;

% Va/pitch (PI)
omega_Vp     = 1;
kp_Va_pit    = (2*zeta*omega_Vp - a_V1) / (-a_V3);
ki_Va_pit    = omega_Vp^2 / (-a_V3);

fprintf('Gains computed:\n');
fprintf('  Roll:     kp=%.4f  kd=%.4f\n', kp_roll, kd_roll);
fprintf('  Course:   kp=%.4f  ki=%.4f\n', kp_course, ki_course);
fprintf('  Pitch:    kp=%.4f  kd=%.4f  KthDC=%.4f\n', kp_pitch, kd_pitch, K_theta_DC);
fprintf('  Altitude: kp=%.4f  ki=%.4f\n', kp_alt, ki_alt);
fprintf('  Va/Thr:   kp=%.4f  ki=%.4f\n', kp_Va_thr, ki_Va_thr);
fprintf('  Va/Pit:   kp=%.4f  ki=%.4f\n', kp_Va_pit, ki_Va_pit);

%% ========================================================================
%%  SIMULATION SETUP
%% ========================================================================
Ts   = SIM.ts_simulation;   % 0.02 s
Tend = 60;                   % simulation duration (s)
t    = 0:Ts:Tend;
N    = length(t);

disturbance_mag  = 0.2;      % rad step disturbance
disturbance_time = 15;       % seconds when disturbance hits

%% ========================================================================
%%  LATERAL AUTOPILOT SIMULATION
%% ========================================================================
% Transfer function model (cascade):
%   phi(s)/delta_a(s) = a_phi2 / (s^2 + a_phi1*s)    [roll plant]
%   chi(s)/phi(s)     = (g/Va)/s                       [course plant]
%   beta(s)/delta_r(s)= a_beta2/(s + a_beta1)          [sideslip plant]

% Discretize plant transfer functions
s_tf = tf('s');
Gphi  = c2d(tf(a_phi2, [1, a_phi1, 0]), Ts, 'zoh');
Gchi  = c2d(tf(g/Va_trim, [1, 0]), Ts, 'zoh');
Gbeta = c2d(tf(a_beta2, [1, a_beta1]), Ts, 'zoh');

% -- Lateral simulation (discrete-time, step commands) --
phi_cmd_val  = 15*pi/180;    % 15 deg step roll command
chi_cmd_val  = 20*pi/180;    % 20 deg step course command
beta_cmd_val = 0;            % zero sideslip commanded

% State arrays
phi_arr       = zeros(1,N);  phi_dot_arr  = zeros(1,N);
chi_arr       = zeros(1,N);
beta_arr      = zeros(1,N);
delta_a_arr   = zeros(1,N);
delta_r_arr   = zeros(1,N);

% Integrators & states
int_chi  = 0;  int_beta = 0;
phi_state  = zeros(2,1);   % [phi; phi_dot] for discrete plant
chi_state  = zeros(1,1);
beta_state = zeros(1,1);

phi_val  = 0; phi_dot_val = 0;
chi_val  = 0; beta_val    = 0;

% Discrete state matrices for phi plant (s^2 + a_phi1*s = s*(s+a_phi1))
[phi_ss_A, phi_ss_B, phi_ss_C, phi_ss_D] = ssdata(Gphi);
phi_x = zeros(size(phi_ss_A,1),1);

[chi_ss_A, chi_ss_B, chi_ss_C, chi_ss_D] = ssdata(Gchi);
chi_x = zeros(size(chi_ss_A,1),1);

[beta_ss_A, beta_ss_B, beta_ss_C, beta_ss_D] = ssdata(Gbeta);
beta_x = zeros(size(beta_ss_A,1),1);

% Simulate lateral loops
phi_cmd_NO_dist = zeros(1,N);     % no disturbance
phi_arr_NO = zeros(1,N);
delta_a_NO = zeros(1,N);
phi_x_nd   = zeros(size(phi_ss_A,1),1);
int_chi_nd = 0;
chi_x_nd   = zeros(size(chi_ss_A,1),1);
chi_arr_nd = zeros(1,N);
phi_cmd_nd = zeros(1,N);
delta_a_arr_nd = zeros(1,N);

for k = 1:N
    tk = t(k);
    
    % --- Step commands (step at t=1s) ---
    if tk >= 1.0
        chi_cmd = chi_cmd_val;
    else
        chi_cmd = 0;
    end
    
    % --- Course → Roll command (PI) ---
    e_chi = chi_cmd - chi_arr_nd(max(1,k-1));
    int_chi_nd = int_chi_nd + Ts*e_chi;
    phi_cmd_nd(k) = sat_val(kp_course*e_chi + ki_course*int_chi_nd, phi_max, -phi_max);
    
    % --- Roll PD ---
    phi_now_nd = phi_arr_NO(max(1,k-1));
    p_now_nd   = (k>1)*(phi_arr_NO(k-1)-phi_arr_NO(max(1,k-2)))/Ts;
    e_phi_nd   = phi_cmd_nd(k) - phi_now_nd;
    da_nd = sat_val(kp_roll*e_phi_nd - kd_roll*p_now_nd, delta_a_max, -delta_a_max);
    delta_a_arr_nd(k) = da_nd;
    
    % --- Roll plant (discrete) ---
    phi_x_nd = phi_ss_A*phi_x_nd + phi_ss_B*da_nd;
    phi_arr_NO(k) = phi_ss_C*phi_x_nd + phi_ss_D*da_nd;
    
    % --- Course plant (discrete) ---
    chi_x_nd = chi_ss_A*chi_x_nd + chi_ss_B*phi_arr_NO(k);
    chi_arr_nd(k) = chi_ss_C*chi_x_nd + chi_ss_D*phi_arr_NO(k);
end

% Simulation WITH disturbance
phi_arr_dist  = zeros(1,N);
chi_arr_dist  = zeros(1,N);
phi_cmd_dist  = zeros(1,N);
delta_a_dist  = zeros(1,N);
int_chi_d = 0;
phi_x_d   = zeros(size(phi_ss_A,1),1);
chi_x_d   = zeros(size(chi_ss_A,1),1);

for k = 1:N
    tk = t(k);
    if tk >= 1.0, chi_cmd = chi_cmd_val; else, chi_cmd = 0; end
    
    e_chi = chi_cmd - chi_arr_dist(max(1,k-1));
    int_chi_d = int_chi_d + Ts*e_chi;
    phi_cmd_dist(k) = sat_val(kp_course*e_chi + ki_course*int_chi_d, phi_max, -phi_max);
    
    phi_now_d = phi_arr_dist(max(1,k-1));
    p_now_d   = (k>1)*(phi_arr_dist(k-1)-phi_arr_dist(max(1,k-2)))/Ts;
    e_phi_d   = phi_cmd_dist(k) - phi_now_d;
    
    % Add disturbance as step at disturbance_time
    dist = 0;
    if tk >= disturbance_time && tk < disturbance_time+Ts
        dist = disturbance_mag;
    end
    
    da_d = sat_val(kp_roll*e_phi_d - kd_roll*p_now_d + dist, delta_a_max, -delta_a_max);
    delta_a_dist(k) = da_d;
    
    phi_x_d = phi_ss_A*phi_x_d + phi_ss_B*da_d;
    phi_arr_dist(k) = phi_ss_C*phi_x_d + phi_ss_D*da_d;
    
    chi_x_d = chi_ss_A*chi_x_d + chi_ss_B*phi_arr_dist(k);
    chi_arr_dist(k) = chi_ss_C*chi_x_d + chi_ss_D*phi_arr_dist(k);
end

% Sideslip simulation (delta_r = -kp_yaw*r, closed loop yaw damper)
[beta_ss_A_cl, beta_ss_B_cl, beta_ss_C_cl, beta_ss_D_cl] = deal([], [], [], []);
beta_arr_sim = zeros(1,N);
delta_r_arr  = zeros(1,N);
int_beta_s   = 0;
beta_x_s     = zeros(size(beta_ss_A,1),1);

for k = 1:N
    tk = t(k);
    if tk >= 1.0, beta_cmd = beta_cmd_val; else, beta_cmd = 0; end
    
    beta_now = beta_arr_sim(max(1,k-1));
    e_beta   = beta_cmd - beta_now;
    int_beta_s = int_beta_s + Ts*e_beta;
    dr = sat_val(kp_sideslip*e_beta + ki_sideslip*int_beta_s, delta_r_max, -delta_r_max);
    delta_r_arr(k) = dr;
    
    beta_x_s = beta_ss_A*beta_x_s + beta_ss_B*dr;
    beta_arr_sim(k) = beta_ss_C*beta_x_s + beta_ss_D*dr;
end

%% ---- LATERAL PLOTS ----------------------------------------------------

figure('Name','Q5 Lateral Autopilot — No Disturbance','NumberTitle','off');
set(gcf,'Position',[50 600 1200 600]);

subplot(2,3,1)
plot(t, phi_arr_NO*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, phi_cmd_nd*180/pi, 'r--', 'LineWidth',1.5);
legend('\phi actual','\phi commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\phi (deg)','FontSize',11);
title('Roll Angle — No Disturbance','FontSize',12); grid on;

subplot(2,3,2)
plot(t, chi_arr_nd*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, repmat(chi_cmd_val*180/pi,1,N).*(t>=1), 'r--', 'LineWidth',1.5);
legend('\chi actual','\chi commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\chi (deg)','FontSize',11);
title('Course Angle','FontSize',12); grid on;

subplot(2,3,3)
plot(t, beta_arr_sim*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, zeros(1,N), 'r--', 'LineWidth',1.5);
legend('\beta actual','\beta commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\beta (deg)','FontSize',11);
title('Sideslip Angle','FontSize',12); grid on;

subplot(2,3,4)
plot(t, delta_a_arr_nd*180/pi, 'b-', 'LineWidth',1.5);
xlabel('Time (s)','FontSize',11); ylabel('\delta_a (deg)','FontSize',11);
title('Aileron Deflection','FontSize',12); grid on;
yline(delta_a_max*180/pi,'r--'); yline(-delta_a_max*180/pi,'r--');

subplot(2,3,5)
plot(t, delta_r_arr*180/pi, 'b-', 'LineWidth',1.5);
xlabel('Time (s)','FontSize',11); ylabel('\delta_r (deg)','FontSize',11);
title('Rudder Deflection','FontSize',12); grid on;
yline(delta_r_max*180/pi,'r--'); yline(-delta_r_max*180/pi,'r--');

sgtitle('Q5: Lateral Autopilot — No Disturbance','FontSize',14,'FontWeight','bold');

figure('Name','Q5 Lateral Autopilot — With Disturbance','NumberTitle','off');
set(gcf,'Position',[50 50 1200 500]);

subplot(1,3,1)
plot(t, phi_arr_dist*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, phi_cmd_dist*180/pi, 'r--', 'LineWidth',1.5);
xline(disturbance_time,'k:','Disturbance','FontSize',9);
legend('\phi actual','\phi commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\phi (deg)','FontSize',11);
title('Roll Angle — With Disturbance (step = 0.2 rad)','FontSize',12); grid on;

subplot(1,3,2)
plot(t, chi_arr_dist*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, repmat(chi_cmd_val*180/pi,1,N).*(t>=1), 'r--', 'LineWidth',1.5);
xline(disturbance_time,'k:','FontSize',9);
legend('\chi actual','\chi commanded','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\chi (deg)','FontSize',11);
title('Course Angle — With Disturbance','FontSize',12); grid on;

subplot(1,3,3)
plot(t, delta_a_dist*180/pi, 'b-', 'LineWidth',1.5);
xlabel('Time (s)','FontSize',11); ylabel('\delta_a (deg)','FontSize',11);
title('Aileron — With Disturbance','FontSize',12); grid on;
yline(delta_a_max*180/pi,'r--'); yline(-delta_a_max*180/pi,'r--');

sgtitle('Q5: Lateral Autopilot — With Disturbance','FontSize',14,'FontWeight','bold');

%% ========================================================================
%%  LONGITUDINAL AUTOPILOT SIMULATION
%% ========================================================================
% Plants:
%   theta(s)/delta_e(s) = a_theta3 / (s^2 + a_theta1*s + a_theta2)
%   h(s)/theta(s)       = Va/s
%   Va(s)/delta_t(s)    = a_V2/(s + a_V1)
%   Va(s)/theta(s)      = -a_V3/(s + a_V1)

Gtheta  = c2d(tf(a_theta3, [1, a_theta1, a_theta2]), Ts, 'zoh');
Gh_th   = c2d(tf(Va_trim, [1, 0]), Ts, 'zoh');
GVa_dt  = c2d(tf(a_V2, [1, a_V1]), Ts, 'zoh');
GVa_th  = c2d(tf(-a_V3, [1, a_V1]), Ts, 'zoh');

[th_A, th_B, th_C, th_D] = ssdata(Gtheta);
[h_A,  h_B,  h_C,  h_D ] = ssdata(Gh_th);
[Vt_A, Vt_B, Vt_C, Vt_D] = ssdata(GVa_dt);
[Vp_A, Vp_B, Vp_C, Vp_D] = ssdata(GVa_th);

theta_cmd_val = 5*pi/180;   % 5 deg pitch step
h_cmd_val     = 50;          % 50 m altitude step
Va_cmd_val    = Va_trim + 10; % 10 m/s airspeed step

theta_arr  = zeros(1,N);
h_arr      = zeros(1,N);
Va_arr_thr = zeros(1,N);   % Va controlled via throttle
Va_arr_pit = zeros(1,N);   % Va controlled via pitch
theta_arr2 = zeros(1,N);   % pitch for Va/pitch loop
delta_e_arr   = zeros(1,N);
delta_t_arr   = zeros(1,N);
delta_e_arr2  = zeros(1,N);
theta_c_arr   = zeros(1,N);
h_arr2        = zeros(1,N);

th_x   = zeros(size(th_A,1),1);
h_x    = zeros(size(h_A,1),1);
Vt_x   = zeros(size(Vt_A,1),1);
Vp_x   = zeros(size(Vp_A,1),1);
th_x2  = zeros(size(th_A,1),1);
h_x2   = zeros(size(h_A,1),1);

int_h  = 0; int_Vt = 0; int_Vp = 0;

% ---- Sim 1: Pitch command tracking ---
th_x_p = zeros(size(th_A,1),1);
int_th  = 0;
theta_cmd_arr = zeros(1,N);

for k = 1:N
    tk = t(k);
    if tk >= 1.0, theta_cmd = theta_cmd_val; else, theta_cmd = 0; end
    theta_cmd_arr(k) = theta_cmd;
    
    theta_now = theta_arr(max(1,k-1));
    q_now = (k>1)*(theta_arr(k-1)-theta_arr(max(1,k-2)))/Ts;
    e_th  = theta_cmd - theta_now;
    de    = sat_val(kp_pitch*e_th - kd_pitch*q_now, delta_e_max, -delta_e_max);
    delta_e_arr(k) = de;
    
    th_x = th_A*th_x + th_B*de;
    theta_arr(k) = th_C*th_x + th_D*de;
end

% ---- Sim 2: Altitude hold (throttle controls Va, pitch controls h) ---
for k = 1:N
    tk = t(k);
    if tk >= 1.0, h_cmd = h_cmd_val; else, h_cmd = 0; end
    if tk >= 1.0, Va_cmd_h = Va_trim; else, Va_cmd_h = Va_trim; end
    
    h_now  = h_arr(max(1,k-1));
    Va_now = Va_arr_thr(max(1,k-1));
    
    % Altitude PI → theta_c
    e_h = h_cmd - h_now;
    int_h = int_h + Ts*e_h;
    theta_c = sat_val(kp_alt*e_h + ki_alt*int_h, theta_c_max, -theta_c_max);
    theta_c_arr(k) = theta_c;
    
    % Pitch PD → delta_e
    theta_now = theta_arr2(max(1,k-1));
    q_now = (k>1)*(theta_arr2(k-1)-theta_arr2(max(1,k-2)))/Ts;
    e_th = theta_c - theta_now;
    de2 = sat_val(kp_pitch*e_th - kd_pitch*q_now, delta_e_max, -delta_e_max);
    delta_e_arr2(k) = de2;
    
    % Airspeed/throttle PI → delta_t
    e_Va = Va_cmd_h - Va_now;
    int_Vt = int_Vt + Ts*e_Va;
    dt = sat_val(kp_Va_thr*e_Va + ki_Va_thr*int_Vt, 1, 0);
    delta_t_arr(k) = dt;
    
    % Update plants
    th_x2 = th_A*th_x2 + th_B*de2;
    theta_arr2(k) = th_C*th_x2 + th_D*de2;
    
    h_x = h_A*h_x + h_B*theta_arr2(k);
    h_arr(k) = h_C*h_x + h_D*theta_arr2(k);
    
    Vt_x = Vt_A*Vt_x + Vt_B*dt;
    dVa_thr = Vt_C*Vt_x + Vt_D*dt;
    Va_arr_thr(k) = Va_trim + dVa_thr;   % perturbation + trim
end

% ---- Sim 3: Airspeed via pitch (Va controlled by theta command) ---
int_Vp2 = 0; th_x3 = zeros(size(th_A,1),1);
h_x3 = zeros(size(h_A,1),1);
Vp_x2 = zeros(size(Vp_A,1),1);
Va_arr_pit = zeros(1,N);

for k = 1:N
    tk = t(k);
    if tk >= 1.0, Va_cmd = Va_cmd_val - Va_trim; else, Va_cmd = 0; end  % perturbation
    
    Va_now = Va_arr_pit(max(1,k-1));
    e_Va = Va_cmd - Va_now;
    int_Vp2 = int_Vp2 + Ts*e_Va;
    theta_c_p = sat_val(kp_Va_pit*e_Va + ki_Va_pit*int_Vp2, theta_c_max, -theta_c_max);
    
    theta_now = theta_arr2(max(1,k-1));  % reuse pitch response
    q_now = (k>1)*(theta_arr2(k-1)-theta_arr2(max(1,k-2)))/Ts;
    de3 = sat_val(kp_pitch*(theta_c_p - theta_now) - kd_pitch*q_now, delta_e_max, -delta_e_max);
    
    th_x3 = th_A*th_x3 + th_B*de3;
    theta_p = th_C*th_x3 + th_D*de3;
    
    Vp_x2 = Vp_A*Vp_x2 + Vp_B*theta_p;
    Va_arr_pit(k) = Vp_C*Vp_x2 + Vp_D*theta_p;
end

%% ---- LONGITUDINAL PLOTS ----------------------------------------------

figure('Name','Q5 Longitudinal — Pitch Tracking','NumberTitle','off');
set(gcf,'Position',[100 600 1200 450]);

subplot(1,3,1)
plot(t, theta_arr*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, theta_cmd_arr*180/pi, 'r--', 'LineWidth',1.5);
legend('\theta actual','\theta commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\theta (deg)','FontSize',11);
title('Pitch Angle Tracking','FontSize',12); grid on;

subplot(1,3,2)
plot(t, h_arr, 'b-', 'LineWidth',1.5); hold on;
plot(t, h_cmd_val*(t>=1), 'r--', 'LineWidth',1.5);
legend('h actual','h commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('Altitude h (m)','FontSize',11);
title('Altitude Hold','FontSize',12); grid on;

subplot(1,3,3)
plot(t, delta_e_arr*180/pi, 'b-', 'LineWidth',1.5); hold on;
plot(t, delta_e_arr2*180/pi, 'm--', 'LineWidth',1.5);
yline(delta_e_max*180/pi,'r:'); yline(-delta_e_max*180/pi,'r:');
legend('\delta_e (pitch cmd)','\delta_e (alt hold)','FontSize',10);
xlabel('Time (s)','FontSize',11); ylabel('\delta_e (deg)','FontSize',11);
title('Elevator Deflection','FontSize',12); grid on;

sgtitle('Q5: Longitudinal — Pitch & Altitude','FontSize',14,'FontWeight','bold');

figure('Name','Q5 Longitudinal — Airspeed Control','NumberTitle','off');
set(gcf,'Position',[100 50 1200 450]);

subplot(1,3,1)
plot(t, Va_arr_thr, 'b-', 'LineWidth',1.5); hold on;
plot(t, Va_trim*ones(1,N), 'r--', 'LineWidth',1.5);
legend('V_a actual','V_a commanded','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('V_a (m/s)','FontSize',11);
title('Airspeed via Throttle','FontSize',12); grid on;

subplot(1,3,2)
plot(t, Va_arr_pit, 'b-', 'LineWidth',1.5); hold on;
plot(t, (Va_cmd_val-Va_trim)*(t>=1), 'r--', 'LineWidth',1.5);
legend('V_a perturbation','V_a commanded (perturb)','FontSize',10,'Location','southeast');
xlabel('Time (s)','FontSize',11); ylabel('\DeltaV_a (m/s)','FontSize',11);
title('Airspeed via Pitch','FontSize',12); grid on;

subplot(1,3,3)
plot(t, delta_t_arr, 'b-', 'LineWidth',1.5);
xlabel('Time (s)','FontSize',11); ylabel('\delta_t','FontSize',11);
title('Throttle Setting','FontSize',12); grid on;
ylim([0 1.1]);

sgtitle('Q5: Longitudinal — Airspeed Control','FontSize',14,'FontWeight','bold');

fprintf('\n=== Q5 Autopilot Simulation Complete ===\n');
fprintf('All plots generated. Check Figure windows.\n');

%% ========================================================================
%% LOCAL UTILITIES
%% ========================================================================
function out = sat_val(in, up, lo)
    out = max(lo, min(up, in));
end

function [x_trim, u_trim, y_trim, climb_rate] = compute_trim_analytic(Va, gamma, R, P)
    g = P.gravity;
    if isinf(R), phi = 0; else, phi = atan(Va^2/(g*abs(R)))*sign(R); end
    if isinf(R), psi_dot = 0; else, psi_dot = Va*cos(gamma)/R; end
    p = 0; q = psi_dot*sin(phi); r = psi_dot*cos(phi);
    n_load = 1/cos(phi);
    W = P.mass*g;
    qbar = 0.5*P.rho*Va^2;
    CL_needed = n_load*W/(qbar*P.S_wing);
    alpha_est = (CL_needed - P.C_L_0)/P.C_L_alpha;
    theta = alpha_est + gamma;
    u_b = Va*cos(alpha_est); v_b = 0; w_b = Va*sin(alpha_est);
    delta_e = -(P.C_m_0 + P.C_m_alpha*alpha_est)/P.C_m_delta_e;
    Drag = qbar*P.S_wing*(P.C_D_0 + P.C_D_alpha*alpha_est + P.C_D_delta_e*delta_e);
    Lift = qbar*P.S_wing*(P.C_L_0 + P.C_L_alpha*alpha_est + P.C_L_delta_e*delta_e);
    T_needed = Drag*cos(alpha_est) - Lift*sin(alpha_est) + W*sin(theta);
    delta_t = max(0, min(1, T_needed/P.T_max));
    beta = 0;
    delta_r = -(P.C_n_0 + P.C_n_beta*beta + P.C_n_r*(P.b/(2*Va))*r)/P.C_n_delta_r;
    delta_a = -(P.C_ell_0 + P.C_ell_beta*beta + P.C_ell_p*(P.b/(2*Va))*p + ...
                P.C_ell_r*(P.b/(2*Va))*r)/P.C_ell_delta_a;
    x_trim   = [0;0;-200; u_b;v_b;w_b; phi;theta;0; p;q;r];
    u_trim   = [delta_e; delta_a; delta_r; delta_t];
    y_trim   = [alpha_est; beta; Va];
    climb_rate = Va*sin(gamma);
end
