function [pt_Z] = calculate_Zpositions_tyres(pc, r, z_offset)
% Computes the positions of the Z-components of one tyre at each
% timestep
% Finds the position on the road which best fits the tyre position
%
% pc is the position of the center of mass (provided in XZY - coordinates)
%           -> Change the order in C++!!! Y-Z
% r - X components of the car body corner position 
%z_offset is the difference between wheel and body due to the spring
%
% Assumes that the car is moving forwards (positive local X velocity)
% 
    n = size(pc, 2);
    
    pt_Z = zeros(1, n);

    % dir is either +1 or -1 -> gives the search direction
    dir = sign(r);
    
    % the tyre is at the same position as the car center (weird but okay)
    if dir == 0
        pt_Z = pc(2,:) + z_offset;
        return;
    end
    
    for i =  1:n
       j = i;
       dist = 0;

       % compare the square of the norm, find the index j, such that the
       % position on the trajectory is beyond the car
       while(dist < r * r)          
           if (j == 1)&&(dir==-1) % reached the beginning of the trajectory
               break;
           end
           if (j == n)&&(dir==1) % reached the end of the trajectory
               break;
           end

           dist_vec = pc(:,i) - pc(:, j);
           dist = dist_vec(1) * dist_vec(1) + dist_vec(3) * dist_vec(3);
           
           j = j + dir;
       end
       
       % linearly interpolate the position of the wheel
       dist_vec = pc(:,i) - pc(:, j-dir);
       dist_inside = sqrt(dist_vec(1) * dist_vec(1) + dist_vec(3) * dist_vec(3));

       dist_vec = pc(:,i) - pc(:, j);
       dist_outside = sqrt(dist_vec(1) * dist_vec(1) + dist_vec(3) * dist_vec(3));

       pt_Z(i) = (dist_inside * pc(2, j-dir) + dist_outside * pc(2, j)) / (dist_inside + dist_outside) + z_offset;

       
    end
end







