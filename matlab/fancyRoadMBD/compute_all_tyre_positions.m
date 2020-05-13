function [pt1, pt2, pt3, pt4, theta, vc, w] = compute_all_tyre_positions(r1, r2, r3, r4, pc, dt)
% rr, rl, fl, fr
%
% pc: list of all XY coordinates along the trajectory [2,:]
% r[1-4]: position of the tyre in the local frame
% dt: timestep
%
%
% pt[1-4]: list of all positions of the tyres [2,:]
% vc: velocity of the center of mass of the car
% theta: orientation angle of the car in the XY plane 
% angular velocity: 
% velocity magnitude from omega + velocity of the center of gravity


% TO NOTICE: 
% MBD only requires the values of pt1, pt2, pt3, pt4, theta(1), vc(:,1), w(1) 
% -> NO NEED TO ITERATE OVER ALL TIMESTEPS for theta, vc, w
%
% ALE only requires the values of pt1(:,1), pt2(:,2), pt4(:,3), pt(:,4), theta, vc(:,1), w(1)
% -> NO NEED TO ITERATE OVER ALL TIMESTEPS for pt1, pt2, pt3, pt4, vc, w
% vc is still somehow required to get all theta


% get the system size
n = size(pc, 2);
inv_dt = 1 / dt;

% add dummy values at the end of the position vector (open end, probably
% sub-optimal)
pc = [pc, pc(:,end)];

% memory allocation
pt1 = zeros(3, n);
pt2 = zeros(3, n);
pt3 = zeros(3, n);
pt4 = zeros(3, n);

theta = zeros(1, n);

vc = zeros(3, n);

w  = zeros(1, n);

% initial velocities (2nd order scheme)
vc(:,1) = (-1.5 * pc(:,1) + 2 * pc(:,2) - 0.5 * pc(:,3)) * inv_dt;

% initial angle
    norm_vcXY = sqrt(vc(1,1)^2 + vc(3,1)^2);
if vc(3,1) >= 0
    theta(1) = acos(vc(1,1) / norm_vcXY);
else
    theta(1) = -acos(vc(1,1) / norm_vcXY);
end    

% initial tyre positions
    c = cos(theta(1));
    s = sin(theta(1));

    pt1(1,1) = pc(1,1) + r1(1) * c - r1(2) * s;
    pt1(3,1) = pc(3,1) + r1(1) * s + r1(2) * c;

    pt2(1,1) = pc(1,1) + r2(1) * c - r2(2) * s;
    pt2(3,1) = pc(3,1) + r2(1) * s + r2(2) * c;
    
    pt3(1,1) = pc(1,1) + r3(1) * c - r3(2) * s;
    pt3(3,1) = pc(3,1) + r3(1) * s + r3(2) * c;
    
    pt4(1,1) = pc(1,1) + r4(1) * c - r4(2) * s;
    pt4(3,1) = pc(3,1) + r4(1) * s + r4(2) * c;

    

    
% update CG velocity, orientation and tyre positions
for i = 2 : n
    % velocity update
    vc(:,i) = (-0.5 * pc(:,i-1) + 0.5 * pc(:,i+1)) * inv_dt;
    
    % angle update
    norm_vcXY = sqrt(vc(1,i)^2 + vc(3,i)^2);
    if norm_vcXY == 0
        theta(i) = theta(i-1);
    elseif vc(3,i) >= 0
        theta(i) = acos(vc(1,i) / norm_vcXY);
    else
        theta(i) = -acos(vc(1,i) / norm_vcXY);
    end 
    
    % case if there is a change in direction
%     if abs(theta(i) - theta(i-1)) > 1
%         theta(i) = theta(i) + pi;
%     end

    % tyre position update
    c = cos(theta(i));
    s = sin(theta(i));

    pt1(1,i) = pc(1,i) + r1(1) * c - r1(2) * s;
    pt1(3,i) = pc(3,i) + r1(1) * s + r1(2) * c;

    pt2(1,i) = pc(1,i) + r2(1) * c - r2(2) * s;
    pt2(3,i) = pc(3,i) + r2(1) * s + r2(2) * c;
    
    pt3(1,i) = pc(1,i) + r3(1) * c - r3(2) * s;
    pt3(3,i) = pc(3,i) + r3(1) * s + r3(2) * c;
    
    pt4(1,i) = pc(1,i) + r4(1) * c - r4(2) * s;
    pt4(3,i) = pc(3,i) + r4(1) * s + r4(2) * c;
        
end

% remove NaN if first theta is NaN -- think of a better way, it is quite
% clumsy 
if isnan(theta(1))
    theta(1) = theta(2);
        c = cos(theta(1));
        s = sin(theta(1));

        pt1(1,1) = pc(1,1) + r1(1) * c - r1(2) * s;
        pt1(3,1) = pc(3,1) + r1(1) * s + r1(2) * c;

        pt2(1,1) = pc(1,1) + r2(1) * c - r2(2) * s;
        pt2(3,1) = pc(3,1) + r2(1) * s + r2(2) * c;

        pt3(1,1) = pc(1,1) + r3(1) * c - r3(2) * s;
        pt3(3,1) = pc(3,1) + r3(1) * s + r3(2) * c;

        pt4(1,1) = pc(1,1) + r4(1) * c - r4(2) * s;
        pt4(3,1) = pc(3,1) + r4(1) * s + r4(2) * c;
end
% what if there is no motion at all ? whole theta will be NaN. TODO





% Part 2: angular velocity update (can probably removed, only the initial w
% is required in MBD or ALE)

% add dummy values at the end of the position vector (open end, probably
% sub-optimal)
theta = [theta, theta(end)];

% initial angular velocity (2nd order scheme)
w(1) = -(-1.5 * theta(1) + 2 * theta(2) - 0.5 * theta(3)) * inv_dt;

% update tyre velocities and angular velocity
for i = 2 : n
    % angular velocity update
    w(i) = -(-0.5 * theta(i-1) + 0.5 * theta(i+1)) * inv_dt;
    
end
end


