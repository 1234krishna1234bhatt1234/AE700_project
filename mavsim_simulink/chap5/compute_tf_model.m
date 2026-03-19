% compute_tf_model.m
% Computes analytical transfer functions (Beard & McLain, Section 5.4)
% for the jet aircraft linearized about the trim point.
%
% Prerequisites: run aerosonde_parameters.m then compute_trim.m first so that
%   P, Va_trim, alpha_trim, theta_trim, gamma_trim, u_trim are in the workspace.
%
% -------------------------------------------------------------------------

% --- Extract trim conditions ---
delta_e_trim = u_trim(1);   % elevator [rad]
delta_t_trim = u_trim(4);   % throttle [0-1]

qbar_trim = 0.5 * P.rho * Va_trim^2;   % dynamic pressure at trim [Pa]

% =========================================================================
% LONGITUDINAL TRANSFER FUNCTION COEFFICIENTS  (Beard & McLain Eq. 5.22-5.30)
% =========================================================================

% Airspeed  Va_dot = -a_V1*(Va-Va*) + a_V2*(dt-dt*) - a_V3*(theta-theta*)
% For jet:  dT/dVa = 0  (thrust independent of airspeed)
C_D_trim = P.C_D_0 + P.C_D_alpha * alpha_trim + P.C_D_delta_e * delta_e_trim;
a_V1 = (P.rho * Va_trim * P.S_wing / P.mass) * C_D_trim;   % [1/s] drag damping
a_V2 = P.T_max / P.mass;                                    % [m/s^2] throttle gain (jet)
a_V3 = P.gravity * cos(theta_trim - gamma_trim);            % [m/s^2/rad] gravity coupling

% Pitch dynamics  theta_ddot + a_th1*theta_dot + a_th2*theta = a_th3*delta_e
a_theta1 = -P.rho * Va_trim * P.S_wing * P.c^2 * P.C_m_q / (4 * P.Jy);     % [1/s]  damping (>0)
a_theta2 = -qbar_trim * P.S_wing * P.c * P.C_m_alpha / P.Jy;               % [1/s^2] stiffness (>0)
a_theta3 =  qbar_trim * P.S_wing * P.c * P.C_m_delta_e / P.Jy;             % [1/s^2/rad] elevator gain

% =========================================================================
% LATERAL TRANSFER FUNCTION COEFFICIENTS  (Beard & McLain Eq. 5.36-5.43)
% =========================================================================

% Roll  phi_ddot + a_phi1*phi_dot = a_phi2*delta_a
a_phi1 = -P.rho * Va_trim * P.S_wing * P.b^2 * P.C_ell_p / (4 * P.Jx);   % [1/s]  roll damping (>0)
a_phi2 =  qbar_trim * P.S_wing * P.b * P.C_ell_delta_a / P.Jx;            % [1/s^2/rad] aileron gain

% Sideslip  beta_dot = -a_beta1*beta + a_beta2*delta_r
a_beta1 = -P.rho * Va_trim * P.S_wing * P.C_Y_beta / (2 * P.mass);         % [1/s]  sideslip stability (>0)
a_beta2 =  P.rho * Va_trim * P.S_wing * P.C_Y_delta_r / (2 * P.mass);      % [1/s/rad] rudder gain

% =========================================================================
% PRINT COEFFICIENTS
% =========================================================================
fprintf('\n--- Longitudinal TF coefficients ---\n');
fprintf('a_V1     = %10.4f  [1/s]\n',        a_V1);
fprintf('a_V2     = %10.4f  [m/s^2 per unit throttle]\n', a_V2);
fprintf('a_V3     = %10.4f  [m/s^2/rad]\n',  a_V3);
fprintf('a_theta1 = %10.4f  [1/s]\n',         a_theta1);
fprintf('a_theta2 = %10.4f  [1/s^2]\n',       a_theta2);
fprintf('a_theta3 = %10.4f  [1/s^2/rad]\n',   a_theta3);

fprintf('\n--- Lateral TF coefficients ---\n');
fprintf('a_phi1   = %10.4f  [1/s]\n',         a_phi1);
fprintf('a_phi2   = %10.4f  [1/s^2/rad]\n',   a_phi2);
fprintf('a_beta1  = %10.4f  [1/s]\n',         a_beta1);
fprintf('a_beta2  = %10.4f  [1/s/rad]\n',     a_beta2);
% =========================================================================
% DEFINE TRANSFER FUNCTIONS  (Beard & McLain Section 5.4)
% =========================================================================

% Longitudinal channel
T_phi_delta_a   = tf([a_phi2],       [1, a_phi1, 0]);          % phi/delta_a
T_chi_phi       = tf([P.gravity/Va_trim], [1, 0]);              % chi/phi  (coordinated turn)
T_theta_delta_e = tf([a_theta3],     [1, a_theta1, a_theta2]); % theta/delta_e
T_h_theta       = tf([Va_trim],      [1, 0]);                   % h/theta
T_h_Va          = tf([theta_trim],   [1, 0]);                   % h/Va
T_Va_delta_t    = tf([a_V2],         [1, a_V1]);                % Va/delta_t
T_Va_theta      = tf([-a_V3],        [1, a_V1]);                % Va/theta

% Lateral channel
T_v_delta_r     = tf([a_beta2],      [1, a_beta1]);             % beta/delta_r

% Display
fprintf('\nTransfer functions created in workspace:\n');
fprintf('  T_phi_delta_a, T_chi_phi, T_theta_delta_e\n');
fprintf('  T_h_theta, T_h_Va, T_Va_delta_t, T_Va_theta, T_v_delta_r\n');

