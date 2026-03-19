function y = autopilot(uu, AP)
%
% autopilot for mavsim
% 
% Modification History:
%   2/11/2010 - RWB
%   5/14/2010 - RWB
%   11/14/2014 - RWB
%   2/16/2019 - RWB
%   
    % process inputs
    NN = 0;
    pn       = uu(1+NN);  % inertial North position
    pe       = uu(2+NN);  % inertial East position
    h        = uu(3+NN);  % altitude
    Va       = uu(4+NN);  % airspeed
    alpha    = uu(5+NN);  % angle of attack
    beta     = uu(6+NN);  % side slip angle
    phi      = uu(7+NN);  % roll angle
    theta    = uu(8+NN);  % pitch angle
    chi      = uu(9+NN);  % course angle
    p        = uu(10+NN); % body frame roll rate
    q        = uu(11+NN); % body frame pitch rate
    r        = uu(12+NN); % body frame yaw rate
    Vg       = uu(13+NN); % ground speed
    wn       = uu(14+NN); % wind North
    we       = uu(15+NN); % wind East
    psi      = uu(16+NN); % heading
    bx       = uu(17+NN); % x-gyro bias
    by       = uu(18+NN); % y-gyro bias
    bz       = uu(19+NN); % z-gyro bias
    NN = NN+19;
    Va_c     = uu(1+NN);  % commanded airspeed (m/s)
    h_c      = uu(2+NN);  % commanded altitude (m)
    chi_c    = uu(3+NN);  % commanded course (rad)
    
%     % If phi_c_ff is specified in Simulink model, then do the following
%     phi_c_ff = uu(4+NN);  % feedforward roll command (rad)
%     NN = NN+4;
    
    % If no phi_c_ff is included in inputs in Simulink model, then do this
    NN = NN+3;
    phi_c_ff = 0;
    
    t        = uu(1+NN);   % time

    % Reset flag for integrators
    if t==0, flag = 0; else, flag = 1; end
    
    %----------------------------------------------------------
    % lateral autopilot
    chi_ref = wrap(chi_c, chi);
    if t==0
        phi_c   = course_with_roll(chi_ref, chi, flag, AP) + phi_c_ff;
        delta_r = yaw_damper(r, flag, AP);
    else
        phi_c   = course_with_roll(chi_ref, chi, flag, AP) + phi_c_ff;
        delta_r = yaw_damper(r, flag, AP);
    end
    delta_a = roll_with_aileron(phi_c, phi, p, AP);  
    
    %----------------------------------------------------------
    % longitudinal autopilot
    
    h_ref = sat(h_c, h+AP.altitude_zone, h-AP.altitude_zone);
    if t==0
        % Altitude state machine logic
        if h <= AP.altitude_takeoff_zone
            delta_t = 1;
            theta_c = 15*pi/180; 
        elseif h < h_c - AP.altitude_hold_zone
            delta_t = 1;
            theta_c = airspeed_with_pitch(Va_c, Va, flag, AP);
        elseif h < h_c + AP.altitude_hold_zone
            delta_t = airspeed_with_throttle(Va_c, Va, flag, AP);
            theta_c = altitude_with_pitch(h_c, h, flag, AP);
        else
            delta_t = 0;
            theta_c = airspeed_with_pitch(Va_c, Va, flag, AP);
        end
    else
        if h <= AP.altitude_takeoff_zone
            delta_t = 1;
            theta_c = 15*pi/180;
        elseif h < h_c - AP.altitude_hold_zone
            delta_t = 1;
            theta_c = airspeed_with_pitch(Va_c, Va, flag, AP);
        elseif h < h_c + AP.altitude_hold_zone
            delta_t = airspeed_with_throttle(Va_c, Va, flag, AP);
            theta_c = altitude_with_pitch(h_c, h, flag, AP);
        else
            delta_t = 0;
            theta_c = airspeed_with_pitch(Va_c, Va, flag, AP);
        end
    end
    delta_e = pitch_with_elevator(theta_c, theta, q, AP);
    
    % limit range of throttle setting to [0,1]
    delta_t = sat(delta_t, 1, 0);
 
    
    %----------------------------------------------------------
    % create outputs
    
    % control outputs
    delta = [delta_e; delta_a; delta_r; delta_t];
    % commanded (desired) states
    x_command = [...
        0;...                    % pn
        0;...                    % pe
        h_c;...                  % h
        Va_c;...                 % Va
        0;...                    % alpha
        0;...                    % beta
        phi_c;...                % phi
        theta_c;...              % theta
        chi_c;...                % chi
        0;...                    % p
        0;...                    % q
        0;...                    % r
        ];
            
    y = [delta; x_command];
 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Autopilot functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phi_c_sat = course_with_roll(chi_c, chi, flag, AP)
    persistent integrator;
    if flag==0, integrator = 0; end
    error = chi_c - chi;
    integrator = integrator + (AP.Ts/2)*error;
    phi_c = AP.course_kp*error + AP.course_ki*integrator;
    phi_c_sat = sat(phi_c, 45*pi/180, -45*pi/180);
end

function delta_a = roll_with_aileron(phi_c, phi, p, AP)
    error = phi_c - phi;
    delta_a = sat(AP.roll_kp*error - AP.roll_kd*p, AP.delta_a_max, -AP.delta_a_max);
end

function delta_e = pitch_with_elevator(theta_c, theta, q, AP)
    error = theta_c - theta;
    delta_e = sat(AP.pitch_kp*error - AP.pitch_kd*q, AP.delta_e_max, -AP.delta_e_max);
end

function delta_t_sat = airspeed_with_throttle(Va_c, Va, flag, AP)
    persistent integrator;
    if flag==0, integrator = 0; end
    error = Va_c - Va;
    integrator = integrator + (AP.Ts/2)*error;
    delta_t = AP.airspeed_throttle_kp*error + AP.airspeed_throttle_ki*integrator;
    delta_t_sat = sat(delta_t, 1, 0);
end

function theta_c_sat = altitude_with_pitch(h_c, h, flag, AP)
    persistent integrator;
    if flag==0, integrator = 0; end
    error = h_c - h;
    integrator = integrator + (AP.Ts/2)*error;
    theta_c = AP.altitude_pitch_kp*error + AP.altitude_pitch_ki*integrator;
    theta_c_sat = sat(theta_c, 30*pi/180, -30*pi/180);
end

function theta_c_sat = airspeed_with_pitch(Va_c, Va, flag, AP)
    persistent integrator;
    if flag==0, integrator = 0; end
    error = Va_c - Va;
    integrator = integrator + (AP.Ts/2)*error;
    theta_c = AP.airspeed_pitch_kp*error + AP.airspeed_pitch_ki*integrator;
    theta_c_sat = sat(theta_c, 30*pi/180, -30*pi/180);
end

function delta_r = yaw_damper(r, flag, AP)
    delta_r = sat(-AP.yaw_damper_kp * r, AP.delta_r_max, -AP.delta_r_max);
end

function out = sat(in, up_limit, low_limit)
    if in > up_limit, out = up_limit;
    elseif in < low_limit, out = low_limit;
    else, out = in; end
end

function chi_c = wrap(chi_c, chi)
    while chi_c - chi > pi, chi_c = chi_c - 2*pi; end
    while chi_c - chi < -pi, chi_c = chi_c + 2*pi; end
end
