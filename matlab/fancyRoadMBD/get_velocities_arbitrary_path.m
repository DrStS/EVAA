function [vc] = get_velocities_arbitrary_path(pc, dt)
% Returns the velocity v at all positions for an arbitrary path described in p
%
% v: all velocities [2, :]
% w: all angular velocities [1, :]
% dt: timestep
%
% v: resultant velocity [2, :]

% get the system size
n = size(pc, 2);
inv_dt = 1 / dt;

% add dummy values at the end of the position vector (open end, probably
% sub-optimal)
pc = [pc, pc(:,end)];

vc = zeros(2, n);

% initial velocities (2nd order scheme)
vc(:,1) = (-1.5 * pc(:,1) + 2 * pc(:,2) - 0.5 * pc(:,3)) * inv_dt;

% update CG velocity, orientation and tyre positions
for i = 2 : n
    % velocity update
    vc(:,i) = (-0.5 * pc(:,i-1) + 0.5 * pc(:,i+1)) * inv_dt;
            
end

end
