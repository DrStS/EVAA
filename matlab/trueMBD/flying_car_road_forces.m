function [FR] = flying_car_road_forces(y, v, m, FT, Froad, d, v_ref, dt)
% y is the current position of the tyre [3]
% v is the current velocity of the tyre [3]
% m is the mass of the tyre
% FT is the force acting on the tyre without any road conditions [3]
% Froad is the force to keep the the tyre on the trajectory [3]
% d is the Z-component of the trajectory [1]
% v_ref is the velocity of the trajectory [3]
% dt is the timestep size
    
    % flying car cases -> no interaction with the road
    FR = FT;
    
    
    
    % if the tyre is on the road (or below), apply upward force
    if y(2)<=d
        FR(1) = Froad(1);
        FR(2) = Froad(2);
        FR(3) = Froad(3);
        if v(2) < 0 % car is moving downwards
            FR(2) = FR(2) - m / dt * v(2);   % !! use higher order method !!
        end
        
        % the road should never pull the car downwards
        if FR(2) < 0
            FR = FT;
        end
        
        % adjust the XY Force accoding to the velocity difference to the trajectory
        % to keep the car on track
        v_diff = v - v_ref;
        v_diff(2) = 0;
        FR = FR - m / dt * v_diff;
        
    end
end

