function [trajectory] = generate_ramp(length_ramp, layout, direction, delta_t, v_init, n)
% always have 10m straight line at the beginning (so that the car is not tilted) 

    direction = direction / norm(direction);
    
    t = linspace(0, delta_t*n, n);

    x = v_init * t;
    dx = v_init * delta_t;

    trajectory(1,:) = x * direction(1);

    trajectory(3,:) = x * direction(2);    

    %ramp between 10 and 10+length
    for i = ceil(10 / dx)+1:floor((10+length_ramp) / dx) 
        trajectory(2,i) = layout(x(i) - 10);
    end
    
end





