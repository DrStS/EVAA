function [] = visualizer3D(fig, y, traj_1, traj_2, plot_traj, delta_t, vel_norms)

set(0,'CurrentFigure',fig)
set(fig,'units','normalized','outerpos',[0 0 1 1.2]);
grid on

num_iter = size(y,1) - 1;

frames_per_second = 5;   % change this to what your graphic card is able
vis_step = ceil(1 / (delta_t * frames_per_second)); % visualization step

% axis handler
% for option 2
axis_update_step = 20 * vis_step;   

% for option 3
scale_factor = 0.8;
minimal_range = 15;
fixed_range_x = minimal_range + scale_factor * vel_norms(1); 
fixed_range_y = minimal_range + scale_factor * vel_norms(1);
fixed_range_z = minimal_range + scale_factor * vel_norms(1);
mid_x = y(1,5);
mid_y = y(1,7);
min_x = mid_x - fixed_range_x / 2;
max_x = mid_x + fixed_range_x / 2;
min_y = mid_y - fixed_range_y / 2;
max_y = mid_y + fixed_range_y / 2;
global_min_z = min(min(min(y(:,33)), min(y(:,36))), min(min(y(:,39)), min(y(:,42)))); 
if plot_traj
    global_min_z = min(min(traj_1(2,:)), min(traj_2(2,:)));
    [road_left, road_right, road_middle, stripe_1, stripe_2, stripe_3, stripe_4, gaps_left, gaps_right] = plotPhysicalRoad(traj_1, traj_2);
    colormap(gray);
    amplitude_gray = max(road_left(3,:)) - min(road_left(3,:));
    axis_gray = [min(road_left(3,:)) - 0.1 * amplitude_gray, max(road_left(3,:)) + 1 * amplitude_gray];
end
axis_min_z = global_min_z;

tic
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
    set(0,'CurrentFigure',fig)
    subplot(1,2,1)
    %car body
    plot3(pcc(1), pcc(2), pcc(3),'gx','LineWidth',2);
     grid on
     hold on
    plot3(pc1(1), pc1(2), pc1(3),'gx','LineWidth',2);
    plot3(pc2(1), pc2(2), pc2(3),'gx','LineWidth',2);
    plot3(pc3(1), pc3(2), pc3(3),'gx','LineWidth',2);
    plot3(pc4(1), pc4(2), pc4(3),'gx','LineWidth',2);
    plot3(  [pc1(1), pcc(1), pc3(1), pc4(1), pcc(1), pc2(1), pc1(1), pc4(1)], ...
            [pc1(2), pcc(2), pc3(2), pc4(2), pcc(2), pc2(2), pc1(2), pc4(2)], ...
            [pc1(3), pcc(3), pc3(3), pc4(3), pcc(3), pc2(3), pc1(3), pc4(3)],'g','LineWidth',2);
    plot3(  [pc2(1), pc3(1)], ...
            [pc2(2), pc3(2)], ...
            [pc2(3), pc3(3)],'g','LineWidth',2);

    %legs
    plot3(pw1(1), pw1(2), pw1(3),'ro','LineWidth',1);
    plot3(pw2(1), pw2(2), pw2(3),'ro','LineWidth',1);
    plot3(pw3(1), pw3(2), pw3(3),'ro','LineWidth',1);
    plot3(pw4(1), pw4(2), pw4(3),'ro','LineWidth',1);

    plot3(pt1(1), pt1(2), pt1(3),'ro','LineWidth',1);
    plot3(pt2(1), pt2(2), pt2(3),'ro','LineWidth',1);
    plot3(pt3(1), pt3(2), pt3(3),'ro','LineWidth',1);
    plot3(pt4(1), pt4(2), pt4(3),'ro','LineWidth',1);

    plot3(  [pc1(1), pw1(1), pt1(1)],...
            [pc1(2), pw1(2), pt1(2)],...
            [pc1(3), pw1(3), pt1(3)],'r','LineWidth',2);

    plot3(  [pc2(1), pw2(1), pt2(1)],...
            [pc2(2), pw2(2), pt2(2)],...
            [pc2(3), pw2(3), pt2(3)],'r','LineWidth',2);

    plot3(  [pc3(1), pw3(1), pt3(1)],...
            [pc3(2), pw3(2), pt3(2)],...
            [pc3(3), pw3(3), pt3(3)],'r','LineWidth',2);

    plot3(  [pc4(1), pw4(1), pt4(1)],...
            [pc4(2), pw4(2), pt4(2)],...
            [pc4(3), pw4(3), pt4(3)],'r','LineWidth',2);

        % road
        if plot_traj
            road_vis_range = max(3,i-10000):min(i+10000,num_iter-3);      
            % road itself
            caxis(axis_gray); 
            s = surf([road_left(1,road_vis_range); stripe_1(1,road_vis_range)],...
                    [road_left(2,road_vis_range); stripe_1(2,road_vis_range)],...
                    [road_left(3,road_vis_range); stripe_1(3,road_vis_range)]);
            s.EdgeColor = 'none';
            
            s = surf([stripe_2(1,road_vis_range); road_middle(1,road_vis_range,1)],...
                    [stripe_2(2,road_vis_range); road_middle(2,road_vis_range,1)],...
                    [stripe_2(3,road_vis_range); road_middle(3,road_vis_range,1)]);
            s.EdgeColor = 'none';
            
            for j=1:size(road_middle,3)-1
                    s = surf([road_middle(1,road_vis_range,j); road_middle(1,road_vis_range,j+1)],...
                             [road_middle(2,road_vis_range,j); road_middle(2,road_vis_range,j+1)],...
                             [road_middle(3,road_vis_range,j); road_middle(3,road_vis_range,j+1)]);
                    s.EdgeColor = 'none';
            end

            s = surf([road_middle(1,road_vis_range,end); stripe_4(1,road_vis_range)],...
                     [road_middle(2,road_vis_range,end); stripe_4(2,road_vis_range)],...
                     [road_middle(3,road_vis_range,end); stripe_4(3,road_vis_range)]);
            s.EdgeColor = 'none';

            
            s = surf([road_right(1,road_vis_range); stripe_3(1,road_vis_range)],...
                    [road_right(2,road_vis_range); stripe_3(2,road_vis_range)],...
                    [road_right(3,road_vis_range); stripe_3(3,road_vis_range)]);
            s.EdgeColor = 'none';
            
        % stripes
        for j = 1:2:length(gaps_right)-2
            range_elems = gaps_right(j):gaps_right(j+1);
            white_right = ones(2,length(range_elems),3);
            
             s = surf([stripe_3(1,range_elems); stripe_4(1,range_elems)],...
                    [stripe_3(2,range_elems); stripe_4(2,range_elems)],...
                    [stripe_3(3,range_elems); stripe_4(3,range_elems)], white_right);
            s.EdgeColor = 'none';
            
            range_elems = gaps_right(j+1):gaps_right(j+2);
            
            if (length(range_elems)>1)
                caxis(axis_gray);
                s = surf([stripe_3(1,range_elems); stripe_4(1,range_elems)],...
                        [stripe_3(2,range_elems); stripe_4(2,range_elems)],...
                        [stripe_3(3,range_elems); stripe_4(3,range_elems)]);
                s.EdgeColor = 'none';
            end

        end
            
        for j = 1:2:length(gaps_left)-2
            range_elems = gaps_left(j):gaps_left(j+1);
            white_left = ones(2,length(range_elems),3);
            
             s = surf([stripe_1(1,range_elems); stripe_2(1,range_elems)],...
                    [stripe_1(2,range_elems); stripe_2(2,range_elems)],...
                    [stripe_1(3,range_elems); stripe_2(3,range_elems)], white_left);
            s.EdgeColor = 'none';
            
            range_elems = gaps_left(j+1):gaps_left(j+2);
            if (length(range_elems)>1)
                caxis(axis_gray);
                s = surf([stripe_1(1,range_elems); stripe_2(1,range_elems)],...
                        [stripe_1(2,range_elems); stripe_2(2,range_elems)],...
                        [stripe_1(3,range_elems); stripe_2(3,range_elems)]);
                s.EdgeColor = 'none';
            end

        end

            
        end
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    timestring=['time=',num2str(floor(delta_t*i/60)),'min',num2str(mod(delta_t*i,60),4),'sec'];
    title(timestring)
     %view(0,0)
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
    min_z = min(min(pt1(3), pt2(3)), min(pt3(3), pt4(3)));

    if ( y(i, 5) + 2.5 > max_x ) 
        % update x-axis
        mid_x = y(i,5);
        fixed_range_x =minimal_range+ scale_factor * vel_norms(i); 
        fixed_range_y =minimal_range+ scale_factor * vel_norms(i);
        fixed_range_z =minimal_range+ scale_factor * vel_norms(i);
        
        min_x = mid_x - 2.5;
        max_x = min_x + fixed_range_x;
        min_y = mid_y - fixed_range_y / 2;
        max_y = mid_y + fixed_range_y / 2;
    end
    
    if ( y(i, 5) - 2.5 < min_x )
        % update x-axis
        mid_x = y(i,5);
        fixed_range_x =minimal_range+ scale_factor * vel_norms(i); 
        fixed_range_y =minimal_range+ scale_factor * vel_norms(i);
        fixed_range_z =minimal_range+ scale_factor * vel_norms(i);
        
        max_x = mid_x + 2.5;
        min_x = max_x - fixed_range_x;
        min_y = mid_y - fixed_range_y / 2;
        max_y = mid_y + fixed_range_y / 2;
    end
    
    if (y(i, 7) + 2.5 > max_y ) 
        % update y-axis
        mid_y = y(i,7);
        fixed_range_x =minimal_range+ scale_factor * vel_norms(i); 
        fixed_range_y =minimal_range+ scale_factor * vel_norms(i);
        fixed_range_z =minimal_range+ scale_factor * vel_norms(i);
        
        min_x = mid_x - fixed_range_x / 2;
        max_x = mid_x + fixed_range_x / 2;
        min_y = mid_y - 2.5;
        max_y = min_y + fixed_range_y;
    end

    if ( y(i, 7) - 2.5 < min_y )
        % update y-axis
        mid_y = y(i,7);
        fixed_range_x =minimal_range+ scale_factor * vel_norms(i); 
        fixed_range_y =minimal_range+ scale_factor * vel_norms(i);
        fixed_range_z =minimal_range+ scale_factor * vel_norms(i);
        
        min_x = mid_x - fixed_range_x / 2;
        max_x = mid_x + fixed_range_x / 2;
        max_y = mid_y + 2.5;
        min_y = max_y - fixed_range_y;
    end
    
    if pcc(3) - fixed_range_z + 1 > axis_min_z
        axis_min_z = pcc(3) - fixed_range_z / 2;
    end
    if pcc(3) - 2 < axis_min_z
        axis_min_z = pcc(3) - fixed_range_z / 2;
    end    
    axis([min_x, max_x, min_y, max_y, axis_min_z-1, axis_min_z + fixed_range_z])

    drawnow;
    
    
    
    
    
    
    set(0,'CurrentFigure',fig)
    subplot(1,2,2)
    
    plot3(pcc(1), pcc(2), pcc(3),'gx','LineWidth',2);
     grid on
     hold on
    plot3(pc1(1), pc1(2), pc1(3),'gx','LineWidth',2);
    plot3(pc2(1), pc2(2), pc2(3),'gx','LineWidth',2);
    plot3(pc3(1), pc3(2), pc3(3),'gx','LineWidth',2);
    plot3(pc4(1), pc4(2), pc4(3),'gx','LineWidth',2);
    plot3(  [pc1(1), pcc(1), pc3(1), pc4(1), pcc(1), pc2(1), pc1(1), pc4(1)], ...
            [pc1(2), pcc(2), pc3(2), pc4(2), pcc(2), pc2(2), pc1(2), pc4(2)], ...
            [pc1(3), pcc(3), pc3(3), pc4(3), pcc(3), pc2(3), pc1(3), pc4(3)],'g','LineWidth',2);
    plot3(  [pc2(1), pc3(1)], ...
            [pc2(2), pc3(2)], ...
            [pc2(3), pc3(3)],'g','LineWidth',2);

    %legs
    plot3(pw1(1), pw1(2), pw1(3),'ro','LineWidth',2);
    plot3(pw2(1), pw2(2), pw2(3),'ro','LineWidth',2);
    plot3(pw3(1), pw3(2), pw3(3),'ro','LineWidth',2);
    plot3(pw4(1), pw4(2), pw4(3),'ro','LineWidth',2);

    plot3(pt1(1), pt1(2), pt1(3),'ro','LineWidth',2);
    plot3(pt2(1), pt2(2), pt2(3),'ro','LineWidth',2);
    plot3(pt3(1), pt3(2), pt3(3),'ro','LineWidth',2);
    plot3(pt4(1), pt4(2), pt4(3),'ro','LineWidth',2);

    plot3(  [pc1(1), pw1(1), pt1(1)],...
            [pc1(2), pw1(2), pt1(2)],...
            [pc1(3), pw1(3), pt1(3)],'r','LineWidth',2);

    plot3(  [pc2(1), pw2(1), pt2(1)],...
            [pc2(2), pw2(2), pt2(2)],...
            [pc2(3), pw2(3), pt2(3)],'r','LineWidth',2);

    plot3(  [pc3(1), pw3(1), pt3(1)],...
            [pc3(2), pw3(2), pt3(2)],...
            [pc3(3), pw3(3), pt3(3)],'r','LineWidth',2);

    plot3(  [pc4(1), pw4(1), pt4(1)],...
            [pc4(2), pw4(2), pt4(2)],...
            [pc4(3), pw4(3), pt4(3)],'r','LineWidth',2);
 
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    velstring = ['speed=',num2str(vel_norms(i)*3.6,3),'km/h'];
    title(velstring)

        % road
    if plot_traj
        road_vis_range = max(1,i-2000):min(i+2000,num_iter);
                
        % road itself
        caxis(axis_gray); 
        s = surf([road_left(1,road_vis_range); stripe_1(1,road_vis_range)],...
                [road_left(2,road_vis_range); stripe_1(2,road_vis_range)],...
                [road_left(3,road_vis_range); stripe_1(3,road_vis_range)]);
        s.EdgeColor = 'none';
            s = surf([stripe_2(1,road_vis_range); road_middle(1,road_vis_range,1)],...
                    [stripe_2(2,road_vis_range); road_middle(2,road_vis_range,1)],...
                    [stripe_2(3,road_vis_range); road_middle(3,road_vis_range,1)]);
            s.EdgeColor = 'none';
            
            for j=1:size(road_middle,3)-1
                    s = surf([road_middle(1,road_vis_range,j); road_middle(1,road_vis_range,j+1)],...
                             [road_middle(2,road_vis_range,j); road_middle(2,road_vis_range,j+1)],...
                             [road_middle(3,road_vis_range,j); road_middle(3,road_vis_range,j+1)]);
                    s.EdgeColor = 'none';
            end

            s = surf([road_middle(1,road_vis_range,end); stripe_4(1,road_vis_range)],...
                     [road_middle(2,road_vis_range,end); stripe_4(2,road_vis_range)],...
                     [road_middle(3,road_vis_range,end); stripe_4(3,road_vis_range)]);
            s.EdgeColor = 'none';


        caxis(axis_gray);
        s = surf([stripe_3(1,road_vis_range); road_right(1,road_vis_range)],...
                [stripe_3(2,road_vis_range); road_right(2,road_vis_range)],...
                [stripe_3(3,road_vis_range); road_right(3,road_vis_range)]);
        s.EdgeColor = 'none';

        % stripes
        for j = 1:2:length(gaps_right)-2
            range_elems = gaps_right(j):gaps_right(j+1);
            white_right = ones(2,length(range_elems),3);
            
             s = surf([stripe_3(1,range_elems); stripe_4(1,range_elems)],...
                    [stripe_3(2,range_elems); stripe_4(2,range_elems)],...
                    [stripe_3(3,range_elems); stripe_4(3,range_elems)], white_right);
            s.EdgeColor = 'none';
            
            range_elems = gaps_right(j+1):gaps_right(j+2);
            
            if (length(range_elems)>1)
                caxis(axis_gray);
                s = surf([stripe_3(1,range_elems); stripe_4(1,range_elems)],...
                        [stripe_3(2,range_elems); stripe_4(2,range_elems)],...
                        [stripe_3(3,range_elems); stripe_4(3,range_elems)]);
                s.EdgeColor = 'none';
            end

        end
            
        for j = 1:2:length(gaps_left)-2
            range_elems = gaps_left(j):gaps_left(j+1);
            white_left = ones(2,length(range_elems),3);
            
             s = surf([stripe_1(1,range_elems); stripe_2(1,range_elems)],...
                    [stripe_1(2,range_elems); stripe_2(2,range_elems)],...
                    [stripe_1(3,range_elems); stripe_2(3,range_elems)], white_left);
            s.EdgeColor = 'none';
            
            range_elems = gaps_left(j+1):gaps_left(j+2);
            if (length(range_elems)>1)
                caxis(axis_gray);
                s = surf([stripe_1(1,range_elems); stripe_2(1,range_elems)],...
                        [stripe_1(2,range_elems); stripe_2(2,range_elems)],...
                        [stripe_1(3,range_elems); stripe_2(3,range_elems)]);
                s.EdgeColor = 'none';
            end

        end
        
        
        plot3(traj_1(1,road_vis_range),...
              traj_1(2,road_vis_range),...
              traj_1(3,road_vis_range), 'b')

          plot3(traj_2(1,road_vis_range),...
                traj_2(2,road_vis_range),...
                traj_2(3,road_vis_range), 'b')
    end
    %view(0,0)
    hold off

    axis([y(i,5) - 3.5, y(i,5) + 3.5, ...
          y(i,7) - 3.5, y(i,7) + 3.5, ...
          min_z-0.5, min_z + 2.5])
    
    drawnow
        
end
vis_time = toc;
showed_frames = floor(num_iter / vis_step);
calculated_frame_rate = showed_frames / vis_time;
vis_string = ['Visualisation took ', num2str(vis_time,4),'sec with a frame rate of ',num2str(calculated_frame_rate,3),' frames/sec.'];
disp(vis_string);
end