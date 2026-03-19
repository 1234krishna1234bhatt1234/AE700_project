% initialize the mav viewer
addpath('../tools');  

% unit conversions
ft_to_m          = 0.3048;
lb_to_kg         = 0.45359237;
slugft2_to_kgm2  = 1.3558179483314;
lb_to_N          = 4.4482216152605;

% initial conditions (set as needed)
MAV.pn0    = 0;    % m
MAV.pe0    = 0;    % m
MAV.pd0    = 0; % m (Down, negative altitude)
MAV.u0     = 200;    % 277 m/s = 1000 km/hr, typical cruise speed 
MAV.v0     = 0;    % m/s
MAV.w0     = 0;    % m/s
MAV.phi0   = 0;    % rad
MAV.theta0 = 0;    % rad
MAV.psi0   = 0;    % rad
MAV.Va0 = 200;   
e = Euler2Quaternion(MAV.phi0, MAV.theta0, MAV.psi0);
MAV.e0     = e(1);
MAV.e1     = e(2);
MAV.e2     = e(3);
MAV.e3     = e(4);
MAV.p0     = 0.2;    % rad/s
MAV.q0     = 0;    % rad/s
MAV.r0     = 0.01;    % rad/s
   
% physical parameters of airframe (from Boeing file)
MAV.gravity = 9.81;
MAV.mass = 636640 * lb_to_kg;
MAV.Jx   = 18200000 * slugft2_to_kgm2;
MAV.Jy   = 33100000 * slugft2_to_kgm2;
MAV.Jz   = 49700000 * slugft2_to_kgm2;
MAV.Jxz  =   970000 * slugft2_to_kgm2;
MAV.S_wing = 5500.0 * (ft_to_m^2);
MAV.b       = 196.0 * ft_to_m;
MAV.c       = 27.3  * ft_to_m;
MAV.rho     = 0.4;      % 0.4 kg/m3 at cruise speed          
MAV.e       = 0.9;          
MAV.AR      = MAV.b^2 / MAV.S_wing;

% Gamma parameters
MAV.Gamma  = MAV.Jx*MAV.Jz - MAV.Jxz^2;
MAV.Gamma1 = (MAV.Jxz*(MAV.Jx-MAV.Jy+MAV.Jz))/MAV.Gamma;
MAV.Gamma2 = (MAV.Jz*(MAV.Jz-MAV.Jy)+MAV.Jxz*MAV.Jxz)/MAV.Gamma;
MAV.Gamma3 = MAV.Jz/MAV.Gamma;
MAV.Gamma4 = MAV.Jxz/MAV.Gamma;
MAV.Gamma5 = (MAV.Jz-MAV.Jx)/MAV.Jy;
MAV.Gamma6 = MAV.Jxz/MAV.Jy;
MAV.Gamma7 = (MAV.Jx*(MAV.Jx-MAV.Jy)+MAV.Jxz*MAV.Jxz)/MAV.Gamma;
MAV.Gamma8 = MAV.Jx/MAV.Gamma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aerodynamic coefficients (per rad, from file; missing terms set to 0)
MAV.C_L_0         = 0.210;
MAV.C_D_0         = 0.0164;
MAV.C_m_0         = 0.0;
MAV.C_L_alpha     = 4.40;
MAV.C_D_alpha     = 0.200;
MAV.C_m_alpha     = -1.0;
MAV.C_L_q         = 6.6;
MAV.C_D_q         = 0.0;
MAV.C_m_q         = -20.5;
MAV.C_L_delta_e   = 0.32;
MAV.C_D_delta_e   = 0;
MAV.C_m_delta_e   = -1.3;
MAV.M             = 0;
MAV.alpha0        = 0;
MAV.epsilon       = 0;
MAV.C_D_p         = 0;

MAV.C_Y_0         = 0;
MAV.C_ell_0       = 0;
MAV.C_n_0         = 0;
MAV.C_Y_beta      = -0.9;
MAV.C_ell_beta    = -0.160;
MAV.C_n_beta      = 0.16;
MAV.C_Y_p         = 0.0;
MAV.C_ell_p       = -0.340;
MAV.C_n_p         = -0.026;
MAV.C_Y_r         = 0;      % not provided
MAV.C_ell_r       = 0.130;
MAV.C_n_r         = -0.280;
MAV.C_Y_delta_a   = 0;
MAV.C_ell_delta_a = -0.013;
MAV.C_n_delta_a   = -0.0018;
MAV.C_Y_delta_r   = 0.120;
MAV.C_ell_delta_r = 0.008;
MAV.C_n_delta_r   = -0.10;

% Engine model (jet, 4× turbofan engines total)
MAV.T_max = 240000 * lb_to_N;   % [N] total max thrust (~1.0755e6 N)
P = MAV;   % alias — all downstream scripts use P
