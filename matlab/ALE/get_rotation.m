function [R] = get_rotation(yaw, pitch, roll)
    R(1,1) = (cos(yaw))*(cos(pitch));
    R(1,2) = ((cos(yaw))*(sin(pitch))*(sin(roll)) - (sin(yaw))*(cos(roll)));
    R(1,3) = ((cos(yaw))*(sin(pitch))*(cos(roll)) + (sin(yaw))*(sin(roll)));
    R(2,1) = (sin(yaw))*(cos(pitch));
    R(2,2) = ((sin(yaw))*(sin(pitch))*(sin(roll)) + (cos(yaw))*(cos(roll)));
    R(2,3) = ((sin(yaw))*(sin(pitch))*(cos(roll)) - (cos(yaw))*(sin(roll)));
    R(3,1) = -(sin(pitch));
    R(3,2) = (cos(pitch))*(sin(roll));
    R(3,3) = (cos(pitch))*(cos(roll));

end