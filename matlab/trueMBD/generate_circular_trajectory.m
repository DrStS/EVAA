function [trajectory] = generate_circular_trajectory(r, position, direction, delta_t, v_init, n)

    direction = direction / norm(direction);
    
    center = position - r * [-direction(2); 0; direction(1)];

    total_distance = v_init * delta_t * n / r;
    
    x = linspace(0, total_distance, n);
    
    trajectory = zeros(3,n);
    
    trajectory(1,:) = center(1) + r * (direction(1) * sin(x) - direction(2) * cos(x));
    trajectory(2,:) = center(2);
    trajectory(3,:) = center(3) + r * (direction(1) * cos(x) + direction(2) * sin(x));

end