% compute trim conditions using 'mavsim_trim.slx'
gamma = 0*pi/180;   % desired flight path angle (radians)
R     = 10000000000;  % desired radius (m): (+) right orbit, (-) left orbit
Va    = 25;

% set initial conditions (12 Euler states: pn,pe,pd,u,v,w,phi,theta,psi,p,q,r)
x0 = [0; 0; -200; Va; 0; 0; 0; 0; 0; 0; 0; 0];
% specify which states to hold equal to the initial conditions
ix = [];

% specify initial inputs
u0 = [...
    0;...  % 1 - delta_e
    0;...  % 2 - delta_a
    0;...  % 3 - delta_r
    1;...  % 4 - delta_t
    ];
% specify which inputs to hold constant
iu = [];

% define constant outputs
y0 = [...
    Va;...  % 1 - Va
    0;...   % 2 - alpha
    0;...   % 3 - beta
    ];
% specify which outputs to hold constant
iy = [1, 3];

% define constant derivatives (12 entries matching 12-state Euler model)
dx0 = [0; 0; -Va*sin(gamma); 0; 0; 0; 0; 0; 0; 0; 0; 0];
%      pn  pe  pd            u  v  w  phi th  psi p  q  r
% pd_dot = -Va*sin(gamma) because pd is positive DOWN, altitude h = -pd
% so climbing (gamma>0) → pd decreasing → pd_dot < 0

if R ~= Inf, dx0(9) = Va*cos(gamma)/R; end  % 9 - psidot
% specify which derivatives to hold constant in trim algorithm
idx = [3; 4; 5; 6; 7; 8; 9; 10; 11; 12];

% compute trim conditions
[x_trim, u_trim, y_trim, dx_trim] = trim('mavsim_trim', x0, u0, y0, ix, iu, iy, dx0, idx);

% check to make sure that the linearization worked (should be small)
norm(dx_trim(3:end) - dx0(3:end))

% set initial conditions to trim conditions
MAV.pn0    = 0;
MAV.pe0    = 0;
MAV.pd0    = -200;
MAV.u0     = x_trim(4);   % u at trim
MAV.v0     = x_trim(5);   % v at trim
MAV.w0     = x_trim(6);   % w at trim
MAV.phi0   = x_trim(7);   % roll angle at trim
MAV.theta0 = x_trim(8);   % pitch angle at trim
MAV.psi0   = x_trim(9);   % yaw angle at trim
MAV.p0     = x_trim(10);  % roll rate at trim
MAV.q0     = x_trim(11);  % pitch rate at trim
MAV.r0     = x_trim(12);  % yaw rate at trim
