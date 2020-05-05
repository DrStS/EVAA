function R = smallAngleApproxMatrix(angles)
    R = zeros(3,3);
    R(1,1) = (cos(angles(2)));
	R(1,2) = ((cos(angles(2))) * (angles(3)) * (angles(1)) - (sin(angles(2))));
	R(1,3) = ((cos(angles(2))) * (angles(3)) + (sin(angles(2))) * (angles(1)));
	R(2,1) = (sin(angles(2)));
	R(2,2) = ((sin(angles(2))) * (angles(3)) * (angles(1)) + (cos(angles(2))));
	R(2,3) = ((sin(angles(2))) * (angles(3)) - (cos(angles(2))) * (angles(1)));
	R(3,1) = -(angles(3));
	R(3,2) = (angles(1));
	R(3,3) = 1;
end