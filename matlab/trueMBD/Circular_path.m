function [FR] = Circular_path(vt, mt, pt)
% vt is the velocity of the tyre
% mt is the mass of the tyre
% pt is the global position of the tyre
%
% The rotation is always around the origin!
    pt(2) = 0;

    radius_tyre = norm(pt);        % curve radius
    
    force_direction_tyre = -pt / radius_tyre;

    velocity_direction_tyre = cross(force_direction_tyre, [0;1;0]);

    velocity_magnitude_tyre = vt'*velocity_direction_tyre;
 
    force_magnitude_tyre = mt * velocity_magnitude_tyre^2 / radius_tyre;
    
    FR = force_magnitude_tyre * force_direction_tyre;
    
end

