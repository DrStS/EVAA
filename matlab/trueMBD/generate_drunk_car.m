function [trajectory] = generate_drunk_car(amplitude, period, delta_t, v_init, n)
   total_distance = v_init * delta_t * n;
    
    x = linspace(0, total_distance, n);
    
    trajectory = zeros(2,n);
    
    trajectory(1,:) = x;
    trajectory(2,:) = amplitude * cos(x * 2 * pi / period);    
    
end
