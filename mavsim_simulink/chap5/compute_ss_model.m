% compute_ss_model.m
% Linearizes the Simulink trim model and extracts the longitudinal and lateral
% decoupled state-space models (Beard & McLain, Equations 5.43 and 5.50).
%
% Prerequisites: run aerosonde_parameters.m, compute_trim.m first.
%
% The 'mavsim_trim' model must have:
%   States : 12-state Euler [pn,pe,pd,u,v,w,phi,theta,psi,p,q,r]
%   Inputs : [delta_e, delta_a, delta_r, delta_t]
%
% -------------------------------------------------------------------------

% 1. Linearize the Simulink model about the trim point
[A, B, C, D] = linmod('mavsim_trim', x_trim, u_trim);

% 2. Extract specific trim states for the kinematic formulas
% States: [pn, pe, pd, u, v, w, phi, theta, psi, p, q, r]
% Indices: 1,  2,  3,  4, 5, 6, 7,   8,     9,   10,11,12
u_star     = x_trim(4);
w_star     = x_trim(6);
phi_star   = x_trim(7);
theta_star = x_trim(8);
p_star     = x_trim(10);
q_star     = x_trim(11);
r_star     = x_trim(12);
g          = 9.81; % Gravity constant

% --- LONGITUDINAL STATE-SPACE (A_lon, B_lon) ---
% Equations for [u_dot; w_dot; q_dot; theta_dot; h_dot]
% Based on your Longitudinal Equation image

A54_lon = u_star*cos(theta_star) + w_star*sin(theta_star);

A_lon = [ A(4,4),  A(4,6),  A(4,11), -g*cos(theta_star), 0;
          A(6,4),  A(6,6),  A(6,11), -g*sin(theta_star), 0;
          A(11,4), A(11,6), A(11,11), 0,                 0;
          0,       0,       1,        0,                 0;
          sin(theta_star), -cos(theta_star), 0, A54_lon, 0 ];

% Inputs: [Elevator (1), Throttle (4)]
% For a jet, throttle acts only in body-x → B(6,4) = B(11,4) = 0 analytically.
% Hardcode those zeros to avoid small linmod numerical artifacts.
B_lon = [ B(4,1),  B(4,4);
          B(6,1),  0;
          B(11,1), 0;
          0,       0;
          0,       0 ];


% --- LATERAL STATE-SPACE (A_lat, B_lat) ---
% Equations for [v_dot; p_dot; r_dot; phi_dot; psi_dot]
% Based on your Lateral Equation image and coefficient table

A14_lat = g * cos(theta_star) * cos(phi_star);
A43_lat = cos(phi_star) * tan(theta_star);
A44_lat = q_star*cos(phi_star)*tan(theta_star) - r_star*sin(phi_star)*tan(theta_star);
A53_lat = cos(phi_star) * sec(theta_star);
A54_lat = p_star*cos(phi_star)*sec(theta_star) - r_star*sin(phi_star)*sec(theta_star);

A_lat = [ A(5,5),  A(5,10), A(5,12), A14_lat, 0;
          A(10,5), A(10,10),A(10,12), 0,       0;
          A(12,5), A(12,10),A(12,12), 0,       0;
          0,       1,       A43_lat, A44_lat, 0;
          0,       0,       A53_lat, A54_lat, 0 ];

% Inputs: [Aileron (2), Rudder (3)]
B_lat = [ B(5,2),  B(5,3);
          B(10,2), B(10,3);
          B(12,2), B(12,3);
          0,       0;
          0,       0 ];