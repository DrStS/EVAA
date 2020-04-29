function [trajectory, numIter] = generate_fancy_road(delta_t)
    
    numIter = 0;
    total_time = 0;
    
    %% accelerated road block (00:00:00)
    acceleration = 0.2;
    v_init = 0.1;
    initial_position = [0; 0; 0];
    initial_direction = [1;0;0];
    time = 1.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    trajectory = generate_accelerated_car(acceleration, initial_position, initial_direction, delta_t, v_init, local_num_iter);
    numIter = numIter + local_num_iter-1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% accelerated road block (00:01.10)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% accelerated road block (00:01.20)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% accelerated road block (00:02.15)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:02.25)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

        %% accelerated road block (00:03.20)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:04.00)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% circular road block (00:04.28)
    radius = -100;
    time = 4.5;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

            %% accelerated road block (00:09.10)
    acceleration = 1;
    time = 1.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
        %% accelerated road block (00:10.20)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% accelerated road block (00:11.00)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% accelerated road block (00:11.25)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:12.05)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

        %% accelerated road block (00:13.10)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:13.20)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% circular road block (00:14.00)
    radius = 100;
    time = 4;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    
        %% accelerated road block (00:18.00)
    acceleration = 0;
    time = 1;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% circular road block (00:19.00)
    radius = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% accelerated road block (00:19.10)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% circular road block (00:20.05)
    radius = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:20.15)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% circular road block (00:21.10)
    radius = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% circular road block (00:21.20)
    radius = 100;
    time = 1.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% circular road block (00:23.15)
    radius = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% accelerated road block (00:23.25)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% circular road block (00:24.20)
    radius = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:25.00)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% circular road block (00:25.25)
    radius = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% accelerated road block (00:26.05)
    acceleration = 1;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% circular road block (00:27.00)
    radius = 10;
    time = 4;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

            %% accelerated road block (00:31.00)
    acceleration = 0.1;
    time = 1;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

        %% circular road block (00:32.00)
    radius = -10;
    time = 4;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

   
            %% accelerated road block (00:36.00)
    acceleration = -0.3;
    time = 1.1;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
        %% jump road block (00:37.05)
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.1 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 1.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
        %% jump road block (00:38.10)
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.1 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 1.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
        %% jump road block (00:39.10)
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.2 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 2.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    

        
    %% jump road block (00:41.20)
        pre_path = 0;
    length_ramp = 5;
    layout = @(X)0.1 * (1 - cos(X * 2 * pi));
    jump = 0;
    time = 2.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

  
                %% accelerated road block (00:43.25)
    acceleration = 0.3;
    time = 2.5;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
     %% jump road block (00:46.00)
    pre_path = 0;
    length_ramp = 4;
    layout = @(X)0.5 * (1 - cos(X/length_ramp * pi));
    jump = 0;
    time = 1.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
   
    
         %% jump road block (00:47.03)
    pre_path = 0;
    length_ramp = 4;
    layout = @(X)0.5 * (cos(X/length_ramp * pi)-1);
    jump = 0;
    time = 1.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


    %% jump road block (00:48.07)
    pre_path = 0;
    length_ramp = 4;
    layout = @(X)0.5 * (1 - cos(X/length_ramp * pi));
    jump = 0;
    time = 1.1;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    
    %%  generate mountain road(00:49.10)    
    acceleration = -0.2;
    layoutZ = @(X)0.5*(cos(X/11 * pi) - 1);
    layoutXY = @(X)2*(1-cos(X/5));
%    layoutXY = @(X)2*(1-sign(cos(X/5)).*sqrt(abs(cos(X/5))));
    time = 3.4;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% accelerated road block (00:52.25)
    acceleration = -2.82;
    layoutZ = @(X)0;
    layoutXY = @(X)-0.05*X.^2;
    time = 1.1;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    direction(2) = 0;

    %% shallow ramp (00:54.00)
    acceleration = 0.8;
    layoutZ = @(X)2*(1-cos(X/30 * pi));
    layoutXY = @(X)-0.001*X.^2;
    time = 8.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
    
    %% turning ascent (1:02.10)
    radius = -11;
    time = 7.6;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
    %% accelerated road block (01:09.25)
    acceleration = -0.2;
    time = 1.2;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% turning ascent (1:10.20)
    radius = 11;
    time = 8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% hills (1:18.25)
    acceleration = 0.6;
    layoutZ = @(X) 0.001 * X.^2 + 1.5 * (cos(X/4.6) - 1);
    layoutXY = @(X)0.0;
    time = 5.2;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


        %% hills (1:24.05)
    acceleration = 0.6;
    layoutZ = @(X) exp(0.005*X)-1 + 2 * (1 - cos(X/10));
    layoutXY = @(X)3 * exp(-X/50) .* (cos(X/12) - 1);
    time = 13;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
       %% straight (1:37.00)
    acceleration = 0;
    layoutZ = @(X) 0;
    layoutXY = @(X) 0;
    time = 0.5;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% jump road block (01:37.20)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X)0;
    jump = -1.5;
    time = 1.46;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    %% jump road block (01:39.05)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X)0;
    jump = -1.5;
    time = 1.46;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    
    
        %% jump road block (01:40.20)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X)0;
    jump = -1.5;
    time = 1.46;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
        %% jump road block (01:42.05)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X)0;
    jump = -1.5;
    time = 1.46;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
        %% slalom (1:43.20)
    acceleration = 0;
    time = 5.86;
    ramp_length = velocity * time + acceleration * time^2;
    layoutZ = @(X) 0.0001* X.^2;
    layoutXY = @(X)(3 * (cos(X/ramp_length * 4 * pi) - 1));
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
%% ascent again (1:49.10)
    acceleration = 2;
    time = 10.5;
    ramp_length = velocity * time + acceleration * time^2;
    layoutZ = @(X) 40 * (X / ramp_length).^2;
    layoutXY = @(X) 10 * (cos(X/ramp_length * pi) - 1);
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    

%% two hills (1:59.25)
    acceleration = 0;
    time = 10.3;
    ramp_length = velocity * time + acceleration * time^2;
    layoutZ = @(X) -10 * (X/ramp_length).^2 + 10 * (1-cos(X/ramp_length * 3 * pi));
    layoutXY = @(X)5 * (1 - cos(X/ramp_length * 4 * pi));
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

   %% jump road block (2:10.05)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:10.20)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9.9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

   %% jump road block (2:11.05)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:11.20)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9.9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
        %% circular road block (02:12.05)
    radius = 200;
    time = 1.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% jump road block (2:14.00)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:14.15)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% jump road block (2:15.00)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:15.15)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
           %% circular road block (02:16.00)
    radius = -200;
    time = 1.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
 
       %% jump road block (2:17.25)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:18.10)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

       %% jump road block (2:18.25)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);


       %% jump road block (2:19.10)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = 9;
    time = 0.5;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    %% last jump (2:19.25)
    acceleration = 5;
    time = 2;
    ramp_length = velocity * time + acceleration * time^2;
    layoutZ = @(X)  30 * (X/ramp_length).^2;
    layoutXY = @(X) 0;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    direction(2) = 0;
    
    %% jump road block (2:21.25)
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X) 0;
    jump = -10;
    time = 2.8;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
   
    direction = direction / norm(direction);
    direction(2) = 0.2;

        %% jump road block (end - (02:24.15))
    pre_path = 0;
    length_ramp = 0.1;
    layout = @(X)0;
    jump = 64.5;
    time = 0.2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);
    
    
    %% end
    time = 2;
    acceleration = -velocity / time;
    ramp_length = velocity * time;
    layoutZ = @(X) 10 * (1-sqrt(cos(X/ramp_length * pi)));
    layoutXY = @(X)0;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_mountain_road(acceleration, layoutXY, layoutZ, position, direction, delta_t, velocity, local_num_iter+1);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = trajectory(:,end)-trajectory(:,end-1);
    position = trajectory(:,end);
    velocity = norm(direction/delta_t);

    
    
end