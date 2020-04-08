function [] = visualizer3D(y, delta_t, vel_norms)
figure()
grid on

num_iter = size(y,1) - 1;

vis_step = 200*max(1, round(1e-3/delta_t)); % visualization step

% axis handler
% for option 2
axis_update_step = 20 * vis_step;   

% for option 3
fixed_range_x = 10 + 0.5 * vel_norms(1); 
fixed_range_y = 10 + 0.5 * vel_norms(1);
fixed_range_z = 1 + 0.05 * vel_norms(1);
mid_x = y(1,5);
mid_y = y(1,7);
min_x = mid_x - fixed_range_x / 2;
max_x = mid_x + fixed_range_x / 2;
min_y = mid_y - fixed_range_y / 2;
max_y = mid_y + fixed_range_y / 2;


for i = 1 : vis_step : num_iter   
    
    %% get components
    pcc = [y(i,5), y(i,7), y(i,6)];
    
    pc1 = [y(i,8), y(i,10), y(i,9)];
    pc2 = [y(i,11), y(i,13), y(i,12)];
    pc3 = [y(i,14), y(i,16), y(i,15)];
    pc4 = [y(i,17), y(i,19), y(i,18)];

    pw1 = [y(i,20), y(i,22), y(i,21)];
    pw2 = [y(i,23), y(i,25), y(i,24)];
    pw3 = [y(i,26), y(i,28), y(i,27)];
    pw4 = [y(i,29), y(i,31), y(i,30)];

    pt1 = [y(i,32), y(i,34), y(i,33)];
    pt2 = [y(i,35), y(i,37), y(i,36)];
    pt3 = [y(i,38), y(i,40), y(i,39)];
    pt4 = [y(i,41), y(i,43), y(i,42)];
 

    
    %% plot
   
    %car body
    plot3(pcc(1), pcc(2), pcc(3),'gx');
     grid on
     hold on
    plot3(pc1(1), pc1(2), pc1(3),'gx');
    plot3(pc2(1), pc2(2), pc2(3),'gx');
    plot3(pc3(1), pc3(2), pc3(3),'gx');
    plot3(pc4(1), pc4(2), pc4(3),'gx');
    plot3(  [pc1(1), pcc(1), pc3(1), pc4(1), pcc(1), pc2(1), pc1(1), pc4(1)], ...
            [pc1(2), pcc(2), pc3(2), pc4(2), pcc(2), pc2(2), pc1(2), pc4(2)], ...
            [pc1(3), pcc(3), pc3(3), pc4(3), pcc(3), pc2(3), pc1(3), pc4(3)],'g');
    plot3(  [pc2(1), pc3(1)], ...
            [pc2(2), pc3(2)], ...
            [pc2(3), pc3(3)],'g');

    %legs
    plot3(pw1(1), pw1(2), pw1(3),'ro');
    plot3(pw2(1), pw2(2), pw2(3),'ro');
    plot3(pw3(1), pw3(2), pw3(3),'ro');
    plot3(pw4(1), pw4(2), pw4(3),'ro');

    plot3(pt1(1), pt1(2), pt1(3),'ro');
    plot3(pt2(1), pt2(2), pt2(3),'ro');
    plot3(pt3(1), pt3(2), pt3(3),'ro');
    plot3(pt4(1), pt4(2), pt4(3),'ro');

    plot3(  [pc1(1), pw1(1), pt1(1)],...
            [pc1(2), pw1(2), pt1(2)],...
            [pc1(3), pw1(3), pt1(3)],'r');

    plot3(  [pc2(1), pw2(1), pt2(1)],...
            [pc2(2), pw2(2), pt2(2)],...
            [pc2(3), pw2(3), pt2(3)],'r');

    plot3(  [pc3(1), pw3(1), pt3(1)],...
            [pc3(2), pw3(2), pt3(2)],...
            [pc3(3), pw3(3), pt3(3)],'r');

    plot3(  [pc4(1), pw4(1), pt4(1)],...
            [pc4(2), pw4(2), pt4(2)],...
            [pc4(3), pw4(3), pt4(3)],'r');
 
    xlabel('X')
    ylabel('Y')
    zlabel('Z')

    
    hold off

    %% axis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% FIRST: fixed proportion, no axis update %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     range_x = max(y(:,5)) - min(y(:,5)) + 5;
%     range_y = max(y(:,7)) - min(y(:,7)) + 5;
%     range = max(range_x, range_y);
%     
%     mid_x = (max(y(:,5)) + min(y(:,5))) / 2;
%     mid_y = (max(y(:,7)) + min(y(:,7))) / 2;
%     
%     min_x = mid_x - range / 2;
%     max_x = mid_x + range / 2;
%     min_y = mid_y - range / 2;
%     max_y = mid_y + range / 2;
%     axis([min_x, max_x, min_y, max_y, -0.4, max(y(:,6))+0.5])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% SECOND: axis update every 20 plots %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     lower_i = max(floor(i / axis_update_step) * axis_update_step, 1);
%     upper_i = min(ceil(i / axis_update_step) * axis_update_step, num_iter);
% 
%     range_x = max(y(lower_i:upper_i, 5)) - min(y(lower_i:upper_i, 5)) + 5;
%     range_y = max(y(lower_i:upper_i, 7)) - min(y(lower_i:upper_i, 7)) + 5;
%     range = max(range_x, range_y);
% 
%     mid_x = (max(y(lower_i:upper_i,5)) + min(y(lower_i:upper_i,5))) / 2;
%     mid_y = (max(y(lower_i:upper_i,7)) + min(y(lower_i:upper_i,7))) / 2;
% 
%     min_x = mid_x - range / 2;
%     max_x = mid_x + range / 2;
%     min_y = mid_y - range / 2;
%     max_y = mid_y + range / 2;
%     
%     axis([min_x, max_x, min_y, max_y, -0.4, max(y(:,6))+0.5])

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% THIRD: fixed proportion, update when necessary %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ( y(i, 5) + 2.5 > max_x ) || ( y(i, 5) - 2.5 < min_x )
        % update x-axis
        mid_x = y(i,5);
        fixed_range_x = 10 + 0.5 * vel_norms(i); 
        fixed_range_y = 10 + 0.5 * vel_norms(i);
        fixed_range_z = 1 + 0.05 * +vel_norms(i);
        
        min_x = mid_x - fixed_range_x / 2;
        max_x = mid_x + fixed_range_x / 2;
        min_y = mid_y - fixed_range_y / 2;
        max_y = mid_y + fixed_range_y / 2;
    end
    
    if (y(i, 7) + 2.5 > max_y ) || ( y(i, 7) - 2.5 < min_y )
        % update x-axis
        mid_y = y(i,7);
        fixed_range_x = 10 + 0.5 * vel_norms(i); 
        fixed_range_y = 10 + 0.5 * vel_norms(i);
        fixed_range_z = 1 + 0.05 * +vel_norms(i);
        
        min_x = mid_x - fixed_range_x / 2;
        max_x = mid_x + fixed_range_x / 2;
        min_y = mid_y - fixed_range_y / 2;
        max_y = mid_y + fixed_range_y / 2;
    end
    
    axis([min_x, max_x, min_y, max_y, -0.4, -0.4 + fixed_range_z])

    drawnow;
        
end
end