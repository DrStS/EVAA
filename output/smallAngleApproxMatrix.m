function R = smallAngleApproxMatrix(angles)
    R = zeros(3,3);
    R(1) = cos(angles(2));
	R(2) = cos(angles(2)) * angles(3) * angles(1) - sin(angles(2));
	R(3) = (cos(angles(2)) * angles(3) + sin(angles(2)) * angles(1));
	R(4) = sin(angles(2));
	R(5) = (sin(angles(2)) * angles(3) * angles(1) + cos(angles(2)));
	R(6) = (sin(angles(2)) * angles(3) - cos(angles(2)) * (angles(1));
	R(7) = -(angles(3));
	R(8) = (angles(1));
	R(9) = 1;
end