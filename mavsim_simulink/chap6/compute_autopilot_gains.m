addpath('../chap5')
load transfer_function_coef
addpath('../parameters')
simulation_parameters

% AP stands for autopilot
AP.gravity = P.gravity;
AP.sigma = 0.05;         % Low-pass filter cutoff for derivative terms
AP.Va0 = P.Va0;          % Nominal airspeed
AP.Ts = P.Ts;            % Sample time

%----------roll loop-------------
% Design for zeta = 0.707 and a fast settling time
zeta_roll = 0.707;
omega_n_roll = 10; % rad/s (Adjust based on aircraft agility)
AP.roll_kp = (omega_n_roll^2) / a_phi2;
AP.roll_kd = (2 * zeta_roll * omega_n_roll - a_phi1) / a_phi2;

%----------course loop-------------
% Course loop is usually 10x slower than roll loop
zeta_course = 0.707;
omega_n_course = omega_n_roll / 10;
AP.course_kp = 2 * zeta_course * omega_n_course * Va / P.gravity;
AP.course_ki = (omega_n_course^2) * Va / P.gravity;

%----------sideslip loop-------------
% Regulates beta to zero using rudder
zeta_beta = 0.707;
omega_n_beta = 5; 
AP.sideslip_kp = (delta_r_max / 20 * pi/180) / (beta_max); % Initial guess
AP.sideslip_ki = 0.1; % Usually very small

%----------yaw damper-------------
% Washes out slow yawing but fights fast oscillations
AP.yaw_damper_tau_r = 0.05; 
AP.yaw_damper_kp = 0.5;    

%----------pitch loop-------------
% Design for zeta = 0.707
zeta_pitch = 0.707;
omega_n_pitch = 10; 
AP.pitch_kp = (omega_n_pitch^2 - a_theta2) / a_theta3;
AP.pitch_kd = (2 * zeta_pitch * omega_n_pitch - a_theta1) / a_theta3;
K_theta_DC = (AP.pitch_kp * a_theta3) / (a_theta2 + AP.pitch_kp * a_theta3);

%----------altitude loop-------------
% Altitude loop is usually 10x slower than pitch loop
zeta_altitude = 0.707;
omega_n_altitude = omega_n_pitch / 10;
AP.altitude_kp = (2 * zeta_altitude * omega_n_altitude) / (K_theta_DC * Va);
AP.altitude_ki = (omega_n_altitude^2) / (K_theta_DC * Va);
AP.altitude_zone = 10;   % Meters: switch between climb/hold/descend

%---------airspeed hold using throttle---------------
zeta_airspeed = 0.707;
omega_n_airspeed = 2; 
AP.airspeed_throttle_kp = (2 * zeta_airspeed * omega_n_airspeed - a_V1) / a_V2;
AP.airspeed_throttle_ki = (omega_n_airspeed^2) / a_V2;