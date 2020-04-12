function [trajectory] = generate_drunk_car(amplitude, period, position, direction, delta_t, v_init, n)

    direction = direction / norm(direction);

    total_distance = v_init * delta_t * n;
    
    x = linspace(0, total_distance, n);
    
    trajectory = zeros(3,n);
    
    path_deviation = amplitude * cos(x * 2 * pi / period);
    
    trajectory(1,:) = position(1) + x * direction(1) - path_deviation * direction(2);
    trajectory(2,:) = position(2);
    trajectory(3,:) = position(3) + x * direction(2) + path_deviation * direction(1);  
end
