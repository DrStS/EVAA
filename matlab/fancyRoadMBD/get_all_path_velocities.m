function [vt1, vt2, vt3, vt4] = get_all_path_velocities(pt1, pt2, pt3, pt4, dt)

    inv_dt = 1 / dt;
    
    n = size(pt1, 2);

    % add dummy values at the end of the position vector (open end, probably
    % sub-optimal)
    pt1 = [pt1, pt1(:,end)];
    pt2 = [pt2, pt2(:,end)];
    pt3 = [pt3, pt3(:,end)];
    pt4 = [pt4, pt4(:,end)];

    
    vt1 = zeros(3, n);
    vt2 = zeros(3, n);
    vt3 = zeros(3, n);
    vt4 = zeros(3, n);


    % initial velocities (2nd order scheme)
    vt1(:,1) = (-1.5 * pt1(:,1) + 2 * pt1(:,2) - 0.5 * pt1(:,3)) * inv_dt;
    vt2(:,1) = (-1.5 * pt2(:,1) + 2 * pt2(:,2) - 0.5 * pt2(:,3)) * inv_dt;
    vt3(:,1) = (-1.5 * pt3(:,1) + 2 * pt3(:,2) - 0.5 * pt3(:,3)) * inv_dt;
    vt4(:,1) = (-1.5 * pt4(:,1) + 2 * pt4(:,2) - 0.5 * pt4(:,3)) * inv_dt;

    for i = 2 : n
    % velocity update
        vt1(:,i) = (-0.5 * pt1(:,i-1) + 0.5 * pt1(:,i+1)) * inv_dt;
        vt2(:,i) = (-0.5 * pt2(:,i-1) + 0.5 * pt2(:,i+1)) * inv_dt;
        vt3(:,i) = (-0.5 * pt3(:,i-1) + 0.5 * pt3(:,i+1)) * inv_dt;
        vt4(:,i) = (-0.5 * pt4(:,i-1) + 0.5 * pt4(:,i+1)) * inv_dt;
    end
end