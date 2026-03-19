function [phi, theta, psi] = Quaternion2Euler(quaternion)
    % converts a quaternion attitude to ZYX Euler angles (roll-pitch-yaw)
    e0 = quaternion(1);
    e1 = quaternion(2);
    e2 = quaternion(3);
    e3 = quaternion(4);

    % roll
    phi = atan2( 2*(e0*e1 + e2*e3), 1 - 2*(e1^2 + e2^2) );
    % pitch
    theta = asin( 2*(e0*e2 - e3*e1) );
    % yaw
    psi = atan2( 2*(e0*e3 + e1*e2), 1 - 2*(e2^2 + e3^2) );
end