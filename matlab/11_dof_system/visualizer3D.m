function [] = visualizer3D(y, delta_t)
figure
grid on
pause(10)

num_iter = size(y,1) - 1;

vis_step = 20*max(1, round(1e-3/delta_t)); % visualization step

for i = 1 : vis_step : num_iter   
    
    %% get components
    pcc = [0, 0 , y(i,5)];

    pc1 = [y(i,18), y(i,22), y(i,6)];
    pc2 = [y(i,19), y(i,23), y(i,7)];
    pc3 = [y(i,20), y(i,24), y(i,8)];
    pc4 = [y(i,21), y(i,25), y(i,9)];

    ps1 = [y(i,18), y(i,22), y(i,10)];
    ps2 = [y(i,19), y(i,23), y(i,11)];
    ps3 = [y(i,20), y(i,24), y(i,12)];
    ps4 = [y(i,21), y(i,25), y(i,13)];

    pw1 = [y(i,18), y(i,22), y(i,14)];
    pw2 = [y(i,19), y(i,23), y(i,15)];
    pw3 = [y(i,20), y(i,24), y(i,16)];
    pw4 = [y(i,21), y(i,25), y(i,17)];
    
    
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
    plot3(ps1(1), ps1(2), ps1(3),'ro');
    plot3(ps2(1), ps2(2), ps2(3),'ro');
    plot3(ps3(1), ps3(2), ps3(3),'ro');
    plot3(ps4(1), ps4(2), ps4(3),'ro');

    plot3(pw1(1), pw1(2), pw1(3),'ro');
    plot3(pw2(1), pw2(2), pw2(3),'ro');
    plot3(pw3(1), pw3(2), pw3(3),'ro');
    plot3(pw4(1), pw4(2), pw4(3),'ro');

    plot3(  [pc1(1), ps1(1), pw1(1)],...
            [pc1(2), ps1(2), pw1(2)],...
            [pc1(3), ps1(3), pw1(3)],'r');

    plot3(  [pc2(1), ps2(1), pw2(1)],...
            [pc2(2), ps2(2), pw2(2)],...
            [pc2(3), ps2(3), pw2(3)],'r');

    plot3(  [pc3(1), ps3(1), pw3(1)],...
            [pc3(2), ps3(2), pw3(2)],...
            [pc3(3), ps3(3), pw3(3)],'r');

    plot3(  [pc4(1), ps4(1), pw4(1)],...
            [pc4(2), ps4(2), pw4(2)],...
            [pc4(3), ps4(3), pw4(3)],'r');
 
    xlabel('X')
    ylabel('Z')
    zlabel('Y')

    axis([min(y(:,18))-0.5 max(y(:,20))+0.5 ...
        min(y(:,23))-0.5 max(y(:,25))+0.5...
        min(min(y(:,14:17)))-0.5 max(max(y(:,5:9)))+0.5])
   
    drawnow;
        
end

end