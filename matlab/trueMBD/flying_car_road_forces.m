function [FR] = flying_car_road_forces(y, v, m, FT, d, dt)
% if the tyre is on the road (or below), apply upward force
    FR = 0;
    if y<=d
        FR = -FT;
        if v < 0 % car is moving downwards
            FR = FR - m * v / dt;   % !! use higher order method !!
        end
    end
end

