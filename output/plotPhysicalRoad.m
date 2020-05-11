function [road_left, road_right, road_middle, stripe_1, stripe_2, stripe_3, stripe_4, gaps_left, gaps_right] = plotPhysicalRoad(traj_left, traj_right)
n = size(traj_left,2);

road_left = zeros(3, n);
road_middle = zeros(3, n, 16);
road_right = zeros(3, n);
stripe_1 = zeros(3, n);
stripe_2 = zeros(3, n);
stripe_3 = zeros(3, n);
stripe_4 = zeros(3, n);

length_right = 0;
length_left = 0;
length_right_white = 0;
length_left_white = 0;

gaps_left = 1;
gaps_right = 1;


gaps = 2;

for i=1:n
    road_left(:,i) = traj_left(:,i) + 0.6 * (traj_left(:,i) - traj_right(:,i));
    road_right(:,i) = traj_right(:,i) + 0.6 * (traj_right(:,i) - traj_left(:,i));
        
    stripe_1(:,i) = traj_left(:,i) + 0.45 * (traj_left(:,i) - traj_right(:,i));
    stripe_2(:,i) = traj_left(:,i) + 0.38 * (traj_left(:,i) - traj_right(:,i));
    
    stripe_3(:,i) = traj_right(:,i) + 0.45 * (traj_right(:,i) - traj_left(:,i));
    stripe_4(:,i) = traj_right(:,i) + 0.38 * (traj_right(:,i) - traj_left(:,i));
    
    for j=1:16
        road_middle(:,i,j) = stripe_2(:,i) - j / 17 * (stripe_2(:,i) - stripe_4(:,i));
    end

    
    % gap handling
    if (i>1)
        diff_right = road_right(:,i) - road_right(:,i-1);
        length_right = length_right + sqrt(diff_right(1) * diff_right(1) + diff_right(2) * diff_right(2));
        length_right_white = length_right_white + sqrt(diff_right(1) * diff_right(1) + diff_right(2) * diff_right(2));

        diff_left = road_left(:,i) - road_left(:,i-1);
        length_left = length_left + sqrt(diff_left(1) * diff_left(1) + diff_left(2) * diff_left(2));
        length_left_white = length_left_white + sqrt(diff_left(1) * diff_left(1) + diff_left(2) * diff_left(2));
    end
    
   if length_right > gaps 
       length_right = 0;
       length_right_white = 0;
       gaps_right = [gaps_right i];
   end    
   if length_left > gaps 
       length_left = 0;
       length_left_white = 0;
       gaps_left = [gaps_left i];
   end    
   
   if (length_right_white > 0.8 * gaps)
       length_right_white = 0;
       gaps_right = [gaps_right i];
   end
   
   if (length_left_white > 0.8 * gaps)
       length_left_white = 0;
       gaps_left = [gaps_left i];
   end
   
end
   gaps_right = [gaps_right n n];
   gaps_left = [gaps_left n n];

end
    