function [sys,x0,str,ts,simStateCompliance] = mav_dynamics(t,x,u,flag,MAV)

switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(MAV);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=mdlDerivatives(t,x,u,MAV);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x);

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end
end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(MAV)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 12;
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 12;
sizes.NumInputs      = 6;
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed

sys = simsizes(sizes);

%
% initialize the initial conditions
%
x0  = [...
    MAV.pn0; MAV.pe0; MAV.pd0;...
    MAV.u0; MAV.v0; MAV.w0;...
    MAV.phi0; MAV.theta0; MAV.psi0;... % Euler angles replaces e0-e3
    MAV.p0; MAV.q0; MAV.r0;...
    ];

%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';
end
% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,uu, MAV)
    % Unpack States
    u = x(4); v = x(5); w = x(6);
    phi = x(7); theta = x(8); psi = x(9); % Euler states
    p = x(10); q = x(11); r = x(12);
    
    % Unpack Forces/Moments
    fx = uu(1); fy = uu(2); fz = uu(3);
    ell = uu(4); m = uu(5); n = uu(6);

    % 1. Position Derivatives (Navigation equations using Euler angles)
    % Standard Rotation Matrix Body -> Level (R_b^v)
    cph = cos(phi);  sph = sin(phi);
    cth = cos(theta); sth = sin(theta);
    cps = cos(psi);  sps = sin(psi);

    pndot = (cth*cps)*u + (sph*sth*cps - cph*sps)*v + (cph*sth*cps + sph*sps)*w;
    pedot = (cth*sps)*u + (sph*sth*sps + cph*cps)*v + (cph*sth*sps - sph*cps)*w;
    pddot = -sth*u + sph*cth*v + cph*cth*w;

    % 2. Velocity Derivatives (Force equations)
    udot = r*v - q*w + fx/MAV.mass;
    vdot = p*w - r*u + fy/MAV.mass;
    wdot = q*u - p*v + fz/MAV.mass;

    % 3. Rotational Kinematics (Euler Rate Equations)
    % This is the part that replaces e0dot...e3dot
    phidot   = p + sin(phi)*tan(theta)*q + cos(phi)*tan(theta)*r;
    thetadot = cos(phi)*q - sin(phi)*r;
    psidot   = (sin(phi)/cos(theta))*q + (cos(phi)/cos(theta))*r;

    % 4. Rotational Dynamics (Moment equations)
    pdot = MAV.Gamma1*p*q - MAV.Gamma2*q*r + MAV.Gamma3*ell + MAV.Gamma4*n;
    qdot = MAV.Gamma5*p*r - MAV.Gamma6*(p^2 - r^2) + (1/MAV.Jy)*m;
    rdot = MAV.Gamma7*p*q - MAV.Gamma1*q*r + MAV.Gamma4*ell + MAV.Gamma8*n;

    sys = [pndot; pedot; pddot; udot; vdot; wdot; phidot; thetadot; psidot; pdot; qdot; rdot];
end
% end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];
end
% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x)
    sys = x;
end

% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)

sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;
end
% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];
end 
% end mdlTerminate
