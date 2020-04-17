function [trajectory, numIter] = generate_fancy_road(delta_t)
    
    numIter = 0;
    total_time = 0;
    
    %% accelerated road block (00:00:00)
    acceleration = 0.2;
    v_init = 0.1;
    initial_position = [0; 0; 0];
    initial_direction = [1;0];
    time = 1.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    trajectory = generate_accelerated_car(acceleration, initial_position, initial_direction, delta_t, v_init, local_num_iter);
    numIter = numIter + local_num_iter-1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    %% accelerated road block (00:01.10)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
    %% accelerated road block (00:01.20)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    %% accelerated road block (00:02.15)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:02.25)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

        %% accelerated road block (00:03.20)
    acceleration = 5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:04.00)
    acceleration = -1;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
    %% circular road block (00:04.28)
    radius = -100;
    time = 4.5;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

            %% accelerated road block (00:09.10)
    acceleration = 1;
    time = 1.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
        %% accelerated road block (00:10.20)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
    %% accelerated road block (00:11.00)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    %% accelerated road block (00:11.25)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:12.05)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

        %% accelerated road block (00:13.10)
    acceleration = -5;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:13.20)
    acceleration = 1.5;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
    %% circular road block (00:14.00)
    radius = 100;
    time = 4;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    
        %% accelerated road block (00:18.00)
    acceleration = 0;
    time = 1;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    %% circular road block (00:10.10)
    radius = 10;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
    %% accelerated road block (00:10.25)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);


    %% circular road block (00:11.15)
    radius = -10;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:12.00)
    acceleration = 0;
    time = 0.9;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    %% circular road block (00:12.20)
    radius = 10;
    time = 0.3;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_circular_trajectory(radius, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

       %% accelerated road block (00:13.00)
    acceleration = 0;
    time = 0.8;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

    
        %% accelerated road block (00:13.25)
    acceleration = -0.1;
    time = 1.2;
    total_time = total_time + time;

    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_accelerated_car(acceleration, position, direction, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;
    
    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);
   
        %% jump road block
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.05 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);
    
    
        %% jump road block
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.1 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);
    
    
        %% jump road block
    pre_path = 0;
    length_ramp = 1;
    layout = @(X)0.2 * (1 - cos(X/length_ramp * 2 * pi));
    jump = 0;
    time = 2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);
    

        
    %% jump road block
    pre_path = 0;
    length_ramp = 2;
    layout = @(X)0.2 * (1 - cos(X/length_ramp * pi));
    jump = 0;
    time = 2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);
    

    %% jump road block
    pre_path = 0;
    length_ramp = 2;
    layout = @(X)0.2 * (cos(X/length_ramp * pi) - 1);
    jump = 0;
    time = 2;
    total_time = total_time + time;
    
    local_num_iter = ceil(time/delta_t);
    new_trajectory = generate_ramp(pre_path, length_ramp, layout, position, direction, jump, delta_t, velocity, local_num_iter+1);
    trajectory = [trajectory(:,1:end-1), new_trajectory];
    numIter = numIter + local_num_iter - 1;

    direction = [trajectory(1,end)-trajectory(1,end-1); trajectory(3,end)-trajectory(3,end-1)];
    position = trajectory(:,end);
    velocity = norm((trajectory(:,end)-trajectory(:,end-1))/delta_t);

end