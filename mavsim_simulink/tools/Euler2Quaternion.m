function quaternion = Euler2Quaternion(phi, theta, psi)
    % Converts ZYX Euler angles (phi: roll, theta: pitch, psi: yaw) to quaternion
    cphi   = cos(phi/2);    sphi   = sin(phi/2);
    ctheta = cos(theta/2);  stheta = sin(theta/2);
    cpsi   = cos(psi/2);    spsi   = sin(psi/2);

    e0 = cphi*ctheta*cpsi + sphi*stheta*spsi;
    e1 = sphi*ctheta*cpsi - cphi*stheta*spsi;
    e2 = cphi*stheta*cpsi + sphi*ctheta*spsi;
    e3 = cphi*ctheta*spsi - sphi*stheta*cpsi;

    quaternion = [e0; e1; e2; e3];
end