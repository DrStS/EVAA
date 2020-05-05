function R = smallAngleApproxMatrix(yaw,pitch,roll)
    R = zeros(3,3);
    R(1,1) = (cos(yaw));
    R(1,2) = ((cos(yaw)) * (pitch) + (sin(yaw)) * (roll));
	R(1,3) = ((cos(yaw)) * (pitch) * (roll) - (sin(yaw)));
	R(2,1) = -(pitch);
    R(2,2) = 1;
	R(2,3) = (roll);
    R(3,1) = (sin(yaw));
    R(3,2) = ((sin(yaw)) * (pitch) - (cos(yaw)) * (roll));
	R(3,3) = ((sin(yaw)) * (pitch) * (roll) + (cos(yaw)));
end

%xyz -> xzy