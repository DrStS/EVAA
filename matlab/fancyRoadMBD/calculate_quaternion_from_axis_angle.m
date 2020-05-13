function q = calculate_quaternion_from_axis_angle(rotation_axis, angle)

% normalize the axis
rotation_axis = rotation_axis / norm(rotation_axis);

q = [rotation_axis(1) * sin(angle / 2); ...
     rotation_axis(2) * sin(angle / 2); ...
     rotation_axis(3) * sin(angle / 2); ...
     cos(angle / 2)];

% normalize the quaternion
    q = q / norm(q);

end