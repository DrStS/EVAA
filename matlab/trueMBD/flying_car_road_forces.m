function [FR] = flying_car_road_forces(y, v, m, FT, d, dt)
% if the tyre is on the road (or below), apply upward force
    FR = [0;0;0];
    if y(2)<=d
        FR(2) = -FT(2);
        if v(2) < 0 % car is moving downwards
            FR(2) = FR(2) - m * v(2) / dt;   % !! use higher order method !!
        end
    end
end

