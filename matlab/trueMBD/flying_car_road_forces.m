function [FR] = flying_car_road_forces(y, v, m, FT, Froad, d, dt)
    
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
    end
end

