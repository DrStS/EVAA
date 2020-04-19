function [trajectory] = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, v_init, n)

    direction = direction / norm(direction);

    t = linspace(0, delta_t*n, n);

    x = v_init * t + acceleration * t.^2 / 2;
    
    trajectory = zeros(3,n);
    
    path_deviation = layoutXY(x);
    
    trajectory(1,:) = position(1) + x * direction(1) + path_deviation * direction(3);
    trajectory(2,:) = position(2) + x * direction(2) + layoutZ(x);
    trajectory(3,:) = position(3) + x * direction(3) - path_deviation * direction(1);  

end
