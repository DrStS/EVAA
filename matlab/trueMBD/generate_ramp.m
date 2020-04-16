function [trajectory] = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, v_init, n)
% always have 10m straight line at the beginning (so that the car is not tilted) 

    direction = direction / norm(direction);
    
    t = linspace(0, delta_t*n, n);

    x = v_init * t;
    dx = v_init * delta_t;

    trajectory(1,:) = position(1) + x * direction(1);
    
    trajectory(2,:) = position(2);

    trajectory(3,:) = position(3) + x * direction(2);    

    %ramp between pre_path and pre_path+length
    for i = ceil(pre_path / dx)+1:floor((pre_path+length_ramp) / dx) 
        trajectory(2,i) = position(2) + layout(x(i) - pre_path);
    end
      
    trajectory(2,i:end) = trajectory(2,i) + jump;


    
end





