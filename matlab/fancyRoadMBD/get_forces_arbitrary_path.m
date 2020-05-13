function [F, T] = get_forces_arbitrary_path(pc, theta, mass, I, dt)
% Returns the force F at all positions for an arbitrary path described in p
%
% pc: all XY-positions [3, :]
% theta: all Z-angles [1, :]
% mass: of the body
% I: moment of inertia in Z direction
% dt: timestep
%
% F: resultant force [3, :]
% T: resultant torque [1, :]


% TO NOTICE
% MBD only needs the forces
% ALE needs the forces and the torque (we might even not need the torque at all) 
% (torque calculation is not checked, some errors might arise from that)

% get the system size
n = size(pc, 2);
inv_dt_squared = 1 / (dt*dt);

% add dummy values at the end of the position vector (open end, probably
% sub-optimal)
pc = [pc, pc(:,end)];

theta = [theta, theta(end)];


% memory allocation
F = zeros(3, n);

T = zeros(1, n);


% initial Force (2nd order scheme)
F(:,1) = (2 * pc(:,1) -5 * pc(:,2) + 4 * pc(:,3) - 1 * pc(:,4)) * inv_dt_squared;

% initial Torque
T(1) = (2 * theta(1) -5 * theta(2) + 4 * theta(3) - 1 * theta(4)) * inv_dt_squared;

for i = 2 : n
    % force update
    F(:,i) = (1 * pc(:,i-1) - 2 * pc(:,i) + 1 * pc(:,i+1)) * inv_dt_squared;

    % torque update
    T(i) = (1 * theta(i-1) - 2 * theta(i) + 1 * theta(i+1)) * inv_dt_squared;
end

% only accelaration was calculated so far
F = mass * F;
T = I * T;


end


