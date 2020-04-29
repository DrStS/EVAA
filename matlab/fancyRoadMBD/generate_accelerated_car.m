function [trajectory] = generate_accelerated_car(acceleration, position, direction, delta_t, v_init, n)
 
    direction = direction / norm(direction);

    t = linspace(0, delta_t*n, n);

    x = v_init * t + acceleration * t.^2 / 2;
    
    trajectory = zeros(3,n);
    
    trajectory(1,:) = x * direction(1) + position(1);
    
    trajectory(2,:) = x * direction(2) + position(2);

    trajectory(3,:) = x * direction(3) + position(3);    
end
