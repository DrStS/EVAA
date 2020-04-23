function [FR] = flying_car_road_forces(m, FT, leg_index, Froad, traj, v_traj, i, dt)
% y is the current position of the tyre [3]
% v is the current velocity of the tyre [3]
% m is the mass of the tyre [1]
% FT is the force acting on the tyre without any road conditions [3]
% Froad is the force to keep the the tyre on the trajectory [3, :]
% traj is the full trajectory of the tyre [3,:]
% v_ref is the velocity of the trajectory [3]
% i is the index of the current timestep [1]
% dt is the timestep size

    global previous_solution_vector previous_force_vector
        
    % trajectory parameters at previous timestep 
   traj_prev = traj(2,i-1);
   vtraj_prev = v_traj(:,i-1);
    
    % get car parameters at previous timestep
    v_prev = previous_solution_vector(16+leg_index*3:16+leg_index*3+2)';
    x_prev = previous_solution_vector(47+leg_index*3:47+leg_index*3+2)';
    f_prev = previous_force_vector(:,leg_index);
    
    % current road forces
    Froad = Froad(:,i);
    FR = Froad;
    
        
    % if the tyre is on the road (or below), apply upward force
    eps = 0;
    if x_prev(2)<=traj_prev+eps
        x_diff = traj_prev - x_prev(2);     % push tyre back to road if it is below the road
        FR(2) = FR(2) + m / (dt) * x_diff;  % NOT PHYSICAL! just in case the velocity correction is not enough
                                            % at high velocities, this might occur
        
        if v_prev(2) < min(0, vtraj_prev(2))-eps % car is moving downwards
            v_diff = min(0, vtraj_prev(2)) - v_prev(2);
            FR(2) = FR(2) + m / dt * v_diff;   % !! use higher order method !!
        end
        % the road should never pull the car downwards
        if Froad(2) < f_prev(2)-eps
            FR(2) = FT(2);
        end
        
    else
        % flying car cases -> no interaction with the road in Z direction
        FR(2) = FT(2);
    end
    
    % adjust the XY Force accoding to the velocity difference to the trajectory
    % to keep the car on track (not really physical, just to have a nice
    % video :P)
    v_diff = v_prev - vtraj_prev;
    v_diff(2) = 0;
    FR = FR - m / dt * v_diff;

end

