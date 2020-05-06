clc; clear; close all
aleSolution_filename = "aleSolution.txt";
aleSolution = [];
parameters = [];
traj_fl = [];
traj_fr = [];
traj_rl = [];
traj_rr = [];
plot_traj = false;

if isfile(aleSolution_filename)
    aleSolution = csvread(aleSolution_filename);
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

% Simulation parameters
delta_t = parameters(1);
numIterations = size(aleSolution, 1);

% car parameters
r1 = [parameters(2); 0; -parameters(6)];
r2 = [parameters(3); 0; -parameters(7)];
r3 = [parameters(4); 0; -parameters(8)];
r4 = [parameters(5); 0; -parameters(9)];

upper_spring_length = [parameters(10); parameters(11); parameters(12); parameters(13)];
lower_spring_length = [parameters(14); parameters(15); parameters(16); parameters(17)];

angles = aleSolution(:,31:33);

% calculate quaternion out of angle representation
%q = calculateQuaternions(angles, numIterations);
%aleSolution(:,31:34) = q;

% rearrange vectors and tyre position
adaptedSolution = calculateCorners(aleSolution, r1, r2, r3, r4);


% calculate velocity norms
vel_norms = zeros(1, numIterations);
for i=1:numIterations
    vel_norms(i) = norm(aleSolution(i, 4:6));
end

% visualize
vis_fig = figure(1);
visualizer3D(vis_fig, adaptedSolution, transpose(traj_rr), transpose(traj_rl), false, delta_t, vel_norms);
