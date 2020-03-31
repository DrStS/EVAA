function [] = visualizer3D(y, delta_t)
figure()
grid on

num_iter = size(y,1) - 1;

vis_step = 200*max(1, round(1e-3/delta_t)); % visualization step

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
    hold off
    
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
    ylabel('Z')
    zlabel('Y')
    
   axis([min(y(:,5))-2.5, max(y(:,5))+2.5,...
       min(y(:,7))-2.5, max(y(:,7))+2.5,...
       min(y(:,6))-2.5, max(y(:,6))+2.5])
   
    drawnow;
        
end
end