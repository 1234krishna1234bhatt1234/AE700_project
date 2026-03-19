% forces_moments.m
%   Computes the forces and moments acting on the airframe. 
%
%   Output is
%       F     - forces
%       M     - moments
%       Va    - airspeed
%       alpha - angle of attack
%       beta  - sideslip angle
%       wind  - wind vector in the inertial frame
%
function out = forces_moments(x, delta, wind, P)
    % relabel the inputs
    pn      = x(1);
    pe      = x(2);
    pd      = x(3);
    u       = x(4);
    v       = x(5);
    w       = x(6);
    phi     = x(7);
    theta   = x(8);
    psi     = x(9);
    p       = x(10);
    q       = x(11);
    r       = x(12);
    delta_e = delta(1);
    delta_a = delta(2);
    delta_r = delta(3);
    delta_t = delta(4);
    w_ns    = wind(1); % steady wind - North
    w_es    = wind(2); % steady wind - East
    w_ds    = wind(3); % steady wind - Down
    u_wg    = wind(4); % gust along body x-axis
    v_wg    = wind(5); % gust along body y-axis    
    w_wg    = wind(6); % gust along body z-axis
    
    % --- 1. COMPUTE WIND DATA ---
    % Pre-compute trig functions for rotation matrix
    cphi = cos(phi); sphi = sin(phi);
    cthe = cos(theta); sthe = sin(theta);
    cpsi = cos(psi); spsi = sin(psi);
    
    % Rotation matrix from NED (World) to Body frame
    R_v_to_b = [
        cthe*cpsi, cthe*spsi, -sthe;
        sphi*sthe*cpsi - cphi*spsi, sphi*sthe*spsi + cphi*cpsi, sphi*cthe;
        cphi*sthe*cpsi + sphi*spsi, cphi*sthe*spsi - sphi*cpsi, cphi*cthe
    ];

    % Transform steady wind into body frame and add gusts
    wind_steady_body = R_v_to_b * [w_ns; w_es; w_ds];
    uw = wind_steady_body(1) + u_wg;
    vw = wind_steady_body(2) + v_wg;
    ww = wind_steady_body(3) + w_wg;
    
    % Calculate total wind in NED frame for the output vector
    wind_total_ned = [w_ns; w_es; w_ds] + R_v_to_b' * [u_wg; v_wg; w_wg];
    w_n = wind_total_ned(1);
    w_e = wind_total_ned(2);
    w_d = wind_total_ned(3);
    
    % --- 2. COMPUTE AIR DATA ---
    % Relative velocity of the air over the wings
    ur = u - uw;
    vr = v - vw;
    wr = w - ww;
    
    % Airspeed (Magnitude of relative velocity)
    Va = sqrt(ur^2 + vr^2 + wr^2);
    
    % Angle of Attack (alpha) and Sideslip (beta)
    if Va == 0
        alpha = 0;
        beta = 0;
    else
        alpha = atan2(wr, ur);
        beta = asin(vr / Va);
    end
    
    % --- 3. COMPUTE EXTERNAL FORCES ---
    % Gravity forces mapped to body frame
    fg_x = -P.mass * P.gravity * sthe;
    fg_y =  P.mass * P.gravity * cthe * sphi;
    fg_z =  P.mass * P.gravity * cthe * cphi;
    
    % Dynamic pressure
    qbar = 0.5 * P.rho * Va^2;
    
    if Va == 0
        fa_x = 0; fa_y = 0; fa_z = 0;
    else
        % Lift and Drag formulas (Longitudinal Aerodynamics)
        Lift = qbar * P.S_wing * (P.C_L_0 + P.C_L_alpha * alpha + P.C_L_q * (P.c / (2*Va)) * q + P.C_L_delta_e * delta_e);
        Drag = qbar * P.S_wing * (P.C_D_0 + P.C_D_alpha * alpha + P.C_D_q * (P.c / (2*Va)) * q + P.C_D_delta_e * delta_e);
        
        % Rotate Lift and Drag from Stability Frame to Body Frame
        fa_x = -Drag * cos(alpha) + Lift * sin(alpha);
        fa_z = -Drag * sin(alpha) - Lift * cos(alpha);
        
        % Lateral Aerodynamics
        fa_y = qbar * P.S_wing * (P.C_Y_0 + P.C_Y_beta * beta + P.C_Y_p * (P.b / (2*Va)) * p + P.C_Y_r * (P.b / (2*Va)) * r + P.C_Y_delta_a * delta_a + P.C_Y_delta_r * delta_r);
    end
    
    % Propulsion Force (Simplified model based on textbook)
    fp_x = P.T_max * delta_t;
    fp_y = 0;
    fp_z = 0;


    
    % Total Forces
    Force(1) = fg_x + fa_x + fp_x;
    Force(2) = fg_y + fa_y + fp_y;
    Force(3) = fg_z + fa_z + fp_z;
    
    % --- 4. COMPUTE EXTERNAL TORQUES ---
    if Va == 0
        Torque(1) = 0; Torque(2) = 0; Torque(3) = 0;
    else
        % Roll Moment (l) - Includes adverse motor torque
        Torque(1) = qbar * P.S_wing * P.b * (P.C_ell_0 + P.C_ell_beta * beta + P.C_ell_p * (P.b / (2*Va)) * p + P.C_ell_r * (P.b / (2*Va)) * r + P.C_ell_delta_a * delta_a + P.C_ell_delta_r * delta_r);
        
        % Pitch Moment (m)
        Torque(2) = qbar * P.S_wing * P.c * (P.C_m_0 + P.C_m_alpha * alpha + P.C_m_q * (P.c / (2*Va)) * q + P.C_m_delta_e * delta_e);
        
        % Yaw Moment (n)
        Torque(3) = qbar * P.S_wing * P.b * (P.C_n_0 + P.C_n_beta * beta + P.C_n_p * (P.b / (2*Va)) * p + P.C_n_r * (P.b / (2*Va)) * r + P.C_n_delta_a * delta_a + P.C_n_delta_r * delta_r);
    end
   
    out = [Force'; Torque'; Va; alpha; beta; w_n; w_e; w_d];
end
