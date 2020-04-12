function [trajectory] = generate_circular_trajectory(r, delta_t, v_init, n)
    total_distance = v_init * delta_t * n / r;
    
    x = linspace(0, total_distance, n);
    
    trajectory = zeros(3,n);
    
    trajectory(1,:) = r * cos(x);
    trajectory(3,:) = r * sin(x);

end