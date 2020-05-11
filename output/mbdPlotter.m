close all
mbdSolution_filename = "mbdSolution.txt";
mdbSolution = [];
parameters = [];
traj_fl = [];
traj_fr = [];
traj_rl = [];
traj_rr = [];
plot_traj = false;

if isfile(mbdSolution_filename)
    mdbSolution = csvread(mbdSolution_filename);
    parameters = csvread('simulationParameters.txt');
else
    disp("Cound not generate any visualization! No output file available!");
    return;
end

if isfile('LegFl.txt') && isfile('LegFr.txt') && isfile('LegRl.txt') && isfile('LegRr.txt')
    traj_fl = csvread('LegFl.txt');
    traj_fr = csvread('LegFr.txt');
    traj_rl = csvread('LegRl.txt');
    traj_rr = csvread('LegRr.txt');
    plot_traj= true;
end

    traj_fl(:,3) = traj_fl(:,3) - traj_fl(2,3) + mdbSolution(2,57);
    traj_fr(:,3) = traj_fr(:,3) - traj_fr(2,3) + mdbSolution(2,60);
    traj_rl(:,3) = traj_rl(:,3) - traj_rl(2,3) + mdbSolution(2,54);
    traj_rr(:,3) = traj_rr(:,3) - traj_rr(2,3) + mdbSolution(2,51);

% Simulation parameters
delta_t = parameters(1);
numIterations = size(mdbSolution, 1);

% car parameters
r1 = [parameters(2); 0; parameters(6)];
r2 = [parameters(3); 0; parameters(7)];
r3 = [parameters(4); 0; parameters(8)];
r4 = [parameters(5); 0; parameters(9)];

upper_spring_length = [parameters(10); parameters(11); parameters(12); parameters(13)];
lower_spring_length = [parameters(14); parameters(15); parameters(16); parameters(17)];

% rearrange vectors and tyre position
adaptedSolution = calculateCornersMBD(mdbSolution, r1, r2, r3, r4);

% calculate velocity norms
vel_norms = zeros(1, numIterations);
for i=1:numIterations
    vel_norms(i) = norm(mdbSolution(i, 4:6));
end

% visualize
vis_fig = figure(1);
visualizer3D(vis_fig, adaptedSolution, traj_rr', traj_rl', plot_traj, false,  delta_t, vel_norms);


