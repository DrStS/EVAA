function [vc, vw1, vw2, vw3, vw4, vt1, vt2, vt3, vt4, omega] = Circular_path_initialization(pcc, pt1, pt2, pt3, pt4, v)

v(2) = 0;    % only consider velocities and displacements in the XZ-plane    
pcc(2) = 0;

radius = norm(pcc);

omega = -cross(pcc, v) / (radius * radius);    % angular velocity of the car

vc = cross(omega, pcc);

vt1 = cross(omega, pt1);
vt2 = cross(omega, pt2);
vt3 = cross(omega, pt3);
vt4 = cross(omega, pt4);

vw1 = vt1;
vw2 = vt2;
vw3 = vt3;
vw4 = vt4;

end
