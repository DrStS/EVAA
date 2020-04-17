function [C] = get_corners(corners, angles)
R = get_rotation(angles(3), angles(1), angles(2));
C = R*corners;
end