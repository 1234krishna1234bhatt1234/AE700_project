addpath('../chap5')
load transfer_function_coef       % loads a_phi1, a_phi2, a_theta1..3, a_V1..3, a_beta1..2, Va_trim, Va, theta_trim
addpath('../parameters')
simulation_parameters              % loads SIM struct (Ts, etc.)
% NOTE: P must be in workspace (aircraft parameters). If not: run aerosonde_parameters then P=MAV;

AP.gravity = P.gravity;
AP.sigma   = 0.05;         % low-pass filter cutoff for derivative terms
AP.Va0     = Va_trim;      % nominal trimmed airspeed
AP.Ts      = SIM.ts_simulation;  % sample time

% Surface deflection limits (problem states 40 deg max)
AP.delta_e_max = 40 * pi/180;   % rad
AP.delta_a_max = 40 * pi/180;   % rad
AP.delta_r_max = 40 * pi/180;   % rad

% Altitude state-machine zones
AP.altitude_takeoff_zone = 10;   % m — below this, full throttle/climb
AP.altitude_hold_zone    = 10;   % m — within this of h_c, hold altitude
AP.altitude_zone         = 10;   % m — legacy variable (same as hold zone)

%----------roll loop (B&M Section 6.1.1)-------------
% Plant: phi(s)/delta_a(s) = a_phi2 / (s^2 + a_phi1*s)
% PD design for desired omega_n and zeta:
%   Char. poly: s^2 + (a_phi1 + a_phi2*kd)*s + a_phi2*kp = s^2 + 2ζω_n*s + ω_n^2
%   → kp = ω_n^2 / a_phi2
%   → kd = (2ζω_n - a_phi1) / a_phi2
zeta_roll    = 0.707;
omega_n_roll = 10;          % rad/s — choose fast enough (< actuator BW)
AP.roll_kp   = (omega_n_roll^2) / a_phi2;
AP.roll_kd   = (2 * zeta_roll * omega_n_roll - a_phi1) / a_phi2;

%----------course loop (B&M Section 6.1.2)-------------
% Plant: chi(s)/phi(s) = g/Va / s  →  PI design, ~10x slower than roll
%   Char. poly: s^2 + (g/Va)*kp*s + (g/Va)*ki = s^2 + 2ζω_n*s + ω_n^2
%   → kp = 2ζω_n * Va/g
%   → ki = ω_n^2 * Va/g
zeta_course    = 0.707;
omega_n_course = omega_n_roll / 10;   % ~10x slower than roll loop
AP.course_kp   = 2 * zeta_course * omega_n_course * Va_trim / P.gravity;
AP.course_ki   = (omega_n_course^2) * Va_trim / P.gravity;

%----------sideslip loop (B&M Section 6.1.3)-------------
% Plant: beta(s)/delta_r(s) = a_beta2 / (s + a_beta1)
% PI design:
%   s^2 + (a_beta1 + a_beta2*kp)*s + a_beta2*ki = s^2 + 2ζω_n*s + ω_n^2
%   → kp = (2ζω_n - a_beta1) / a_beta2
%   → ki = ω_n^2 / a_beta2
zeta_beta       = 0.707;
omega_n_beta    = 5;        % rad/s
AP.sideslip_kp  = (2 * zeta_beta * omega_n_beta - a_beta1) / a_beta2;
AP.sideslip_ki  = (omega_n_beta^2) / a_beta2;

%----------yaw damper (B&M Section 6.1.4)-------------
AP.yaw_damper_tau_r = 0.05;   % washout filter time constant
AP.yaw_damper_kp    = 0.5;

%----------pitch loop (B&M Section 6.2.1)-------------
% Plant: theta(s)/delta_e(s) = a_theta3 / (s^2 + a_theta1*s + a_theta2)
% PD design:
%   Char. poly: s^2 + (a_theta1 + a_theta3*kd)*s + (a_theta2 + a_theta3*kp)
%   → kp = (ω_n^2 - a_theta2) / a_theta3
%   → kd = (2ζω_n - a_theta1) / a_theta3
zeta_pitch    = 0.707;
omega_n_pitch = 10;         % rad/s
AP.pitch_kp   = (omega_n_pitch^2 - a_theta2) / a_theta3;
AP.pitch_kd   = (2 * zeta_pitch * omega_n_pitch - a_theta1) / a_theta3;

% DC gain of closed pitch inner loop (used in altitude design):
%   K_theta_DC = a_theta3*kp / (a_theta2 + a_theta3*kp) = ω_n^2/(ω_n^2) = 1 theoretically
K_theta_DC = (AP.pitch_kp * a_theta3) / (a_theta2 + AP.pitch_kp * a_theta3);

%----------altitude loop (B&M Section 6.2.2)-------------
% Plant: h(s)/theta(s) = Va/s   (through pitch → airspeed → altitude)
% With K_theta_DC, effective plant: h(s)/theta_c(s) = K_theta_DC*Va / s
% PI design:
%   → kp = 2ζω_n / (K_theta_DC * Va)
%   → ki = ω_n^2 / (K_theta_DC * Va)
zeta_altitude    = 0.707;
omega_n_altitude = omega_n_pitch / 10;   % ~10x slower than pitch loop
AP.altitude_pitch_kp = (2 * zeta_altitude * omega_n_altitude) / (K_theta_DC * Va_trim);
AP.altitude_pitch_ki = (omega_n_altitude^2) / (K_theta_DC * Va_trim);

%---------airspeed hold using throttle (B&M Section 6.2.3)---------------
% Plant: Va(s)/delta_t(s) = a_V2 / (s + a_V1)
% PI design:
%   s^2 + (a_V1 + a_V2*kp)*s + a_V2*ki = s^2 + 2ζω_n*s + ω_n^2
%   → kp = (2ζω_n - a_V1) / a_V2
%   → ki = ω_n^2 / a_V2
zeta_Va_throttle    = 0.707;
omega_n_Va_throttle = 2;    % rad/s
AP.airspeed_throttle_kp = (2 * zeta_Va_throttle * omega_n_Va_throttle - a_V1) / a_V2;
AP.airspeed_throttle_ki = (omega_n_Va_throttle^2) / a_V2;

%---------airspeed hold using pitch (B&M Section 6.2.4)------------------
% Plant: Va(s)/theta(s) = -a_V3 / (s + a_V1)   [negative gain: pitch up → Va decreases]
% PI design (note negative plant gain → negative kp for correct sign):
%   → kp = (2ζω_n - a_V1) / (-a_V3)
%   → ki = ω_n^2 / (-a_V3)
zeta_Va_pitch    = 0.707;
omega_n_Va_pitch = 1;       % rad/s — must be slower than pitch loop
AP.airspeed_pitch_kp = (2 * zeta_Va_pitch * omega_n_Va_pitch - a_V1) / (-a_V3);
AP.airspeed_pitch_ki = (omega_n_Va_pitch^2) / (-a_V3);

disp('Autopilot gains computed:')
fprintf('  Roll:      kp=%.4f  kd=%.4f\n', AP.roll_kp, AP.roll_kd)
fprintf('  Course:    kp=%.4f  ki=%.4f\n', AP.course_kp, AP.course_ki)
fprintf('  Pitch:     kp=%.4f  kd=%.4f  K_DC=%.4f\n', AP.pitch_kp, AP.pitch_kd, K_theta_DC)
fprintf('  Altitude:  kp=%.4f  ki=%.4f\n', AP.altitude_pitch_kp, AP.altitude_pitch_ki)
fprintf('  Va/Throttle: kp=%.4f  ki=%.4f\n', AP.airspeed_throttle_kp, AP.airspeed_throttle_ki)
fprintf('  Va/Pitch:    kp=%.4f  ki=%.4f\n', AP.airspeed_pitch_kp, AP.airspeed_pitch_ki)
