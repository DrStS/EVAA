function [q, angle, rotation_axis] = get_quaternion(v1, v2)
% this function computes a possible unit quaternion to obtain v2 from v1
% FOR CPP: split this function in two:
%   -one returns only the quaternion 
%   -one returns only angle+axis (remove all parts that involve)
% return the quaternion, the axis and the angle

v1 = v1 / norm(v1);
v2 = v2 / norm(v2);

rotation_axis = cross(v1, v2);


% compute the angle between v1 and v2 (using the dot product) 
angle = acos(dot(v1,v2));

if (~isreal(angle))
    angle = real(angle);
end


% if the crossproduct between v1 and v2 is almost zero, assume both vectors
% are parallel and chose an arbitrary rotation axis orthogonal to v1 
%
%(no need to do it for the axis/angle representation, since this corresponds
% to an angle of 0° (torque == 0))
%
if norm(rotation_axis) < 1e-12
    if v1(3) ~= 0 
        rotation_axis = [-1; -1; (v1(1) + v1(2)) / v1(3) ];
    elseif v1(2) ~= 0
        rotation_axis = [-1; (v1(1) + v1(3)) / v1(2); -1 ];
    elseif v1(1) ~= 0
        rotation_axis = [(v1(3) + v1(2)) / v1(1); -1; -1 ];
    else
        rotation_axis = [1 0 0];
    end
end

rotation_axis = rotation_axis / norm(rotation_axis);


q = [rotation_axis(1) * sin(angle / 2); ...
     rotation_axis(2) * sin(angle / 2); ...
     rotation_axis(3) * sin(angle / 2); ...
     cos(angle / 2)];

% normalize the quaternion
q = q / norm(q);

%% verify if the rotation is correct (this might be quite expensive an seems never to happen, so can probably be removed)
% get rotation matrix R
R = zeros(3);
R(:,1) = [1 - 2 * (q(2)^2 + q(3)^2); ...
          2 * (q(1) * q(2) + q(3) * q(4)); ...
          2 * (q(1) * q(3) - q(2) * q(4))];

R(:,2) = [2 * (q(1) * q(2) - q(3) * q(4)); ...
          1 - 2 * (q(1)^2 + q(3)^2); ...
          2 * (q(2) * q(3) + q(1) * q(4))];

R(:,3) = [2 * (q(1) * q(3) + q(2) * q(4)); ...
          2 * (q(2) * q(3) - q(1) * q(4)); ...
          1 - 2 * (q(1)^2 + q(2)^2)];

      
% check the vector (this part can probably be removed)
if norm(v2 - R * v1) > 1e-7
    q(1:3) = -q(1:3);
    rotation_axis = -rotation_axis; % do we need it? 
end

end