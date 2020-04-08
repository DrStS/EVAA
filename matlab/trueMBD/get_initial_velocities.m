function [vt1, vt2, vt3, vt4] = get_initial_velocities(pt1, pt2, pt3, pt4, dt)

    inv_dt = 1 / dt;

    % initial velocities (2nd order scheme)
    vt1 = (-1.5 * pt1(:,1) + 2 * pt1(:,2) - 0.5 * pt1(:,3)) * inv_dt;
    vt2 = (-1.5 * pt2(:,1) + 2 * pt2(:,2) - 0.5 * pt2(:,3)) * inv_dt;
    vt3 = (-1.5 * pt3(:,1) + 2 * pt3(:,2) - 0.5 * pt3(:,3)) * inv_dt;
    vt4 = (-1.5 * pt4(:,1) + 2 * pt4(:,2) - 0.5 * pt4(:,3)) * inv_dt;

end