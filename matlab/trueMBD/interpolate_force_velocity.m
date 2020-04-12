function [d, v, F] = interpolate_force_velocity(v_traj, F_traj, x_current, traj, i)
% get force and velocity values for any point between two points on a
% trajectory (evaluate distance to trajectory only with XY coordinates)

if (i>1)
    dist_previous = x_current - traj(:,i-1);
    dist_next = traj(:,i) - x_current;
    
    sum_distance = traj(:,i) - traj(:,i-1);
    
    F = zeros(3,1);
    v = zeros(3,1);
    d = zeros(3,1);
    
    eps = 1e-10;
    for j = 1:3
        if abs(sum_distance(j)) > eps
            if dist_next(j) / sum_distance(j) < 0
                F(j) = F_traj(j,i);
                v(j) = v_traj(j,i);
                d(j) = traj(j,i);
            elseif dist_previous(j) / sum_distance(j) < 0
                F(j) = F_traj(j,i-1);
                v(j) = v_traj(j,i-1);
                d(j) = traj(j,i-1);
            else
                F(j) = (dist_next(j) * F_traj(j,i-1) + dist_previous(j) * F_traj(j,i)) / sum_distance(j);
                v(j) = (dist_next(j) * v_traj(j,i-1) + dist_previous(j) * v_traj(j,i)) / sum_distance(j);
                d(j) = (dist_next(j) * traj(j,i-1) + dist_previous(j) * traj(j,i)) / sum_distance(j);
            end
        else
            F(j) = F_traj(j,i);
            v(j) = v_traj(j,i);
            d(j) = traj(j,i);
        end
    end
    
else
    v = v_traj(:,1);
    F = F_traj(:,1);
    d = traj(:,1);
end
end