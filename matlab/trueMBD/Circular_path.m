function [FR] = Circular_path(vt, vc, mc, mt, FT, pt, pcc)
% v is the velocity of the mass body
% m is a fourth of the total mass of the car
% FT is the current force in the tyre
% r is the position of the tyre in the local car frame with respect to the
% center of mass
% pt is the global position of the body
%
% The rotation is always around the origin!

    pt(2) = 0;
    pcc(2) = 0; 
        
    radius = norm(pt);        % curve radius
    radius_body = norm(pcc);  % curve radius
    
    force_direction = -pt / radius;

    force_direction_body = -pcc / radius_body;

%    velocity_direction = cross(force_direction, [0;1;0]);
    
%    velocity_magnitude = v'*velocity_direction;
     velocity_magnitude_body = norm(vc);
     velocity_magnitude = norm(vt);
    
%    force_direction = cross([0;1;0], velocity_direction);
    
    force_magnitude = mt * velocity_magnitude^2 / radius;
    force_magnitude_body = mc * velocity_magnitude_body^2 / radius_body;
    
    FR = force_magnitude*force_direction + force_magnitude_body*force_direction_body; 

    FR(2) = -FT(2);
    
end

