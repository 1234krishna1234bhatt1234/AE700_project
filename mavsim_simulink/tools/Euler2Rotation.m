function R = Euler2Rotation(phi, theta, psi)
    % rotation is body to inertial (NED) frame, using ZYX convention
    cphi = cos(phi);   sphi = sin(phi);
    cth  = cos(theta); sth  = sin(theta);
    cpsi = cos(psi);   spsi = sin(psi);

    R = [ ...
        cth*cpsi,                 cth*spsi,               -sth;
        sphi*sth*cpsi - cphi*spsi, sphi*sth*spsi + cphi*cpsi, sphi*cth;
        cphi*sth*cpsi + sphi*spsi, cphi*sth*spsi - sphi*cpsi, cphi*cth ...
        ];
end