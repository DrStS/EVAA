%
% Only compare two results if the initial spring lengths is the same
%


mbdSolution_filename = "mbdSolution.txt";
mdbSolution = [];

aleSolution_filename = "aleSolution.txt";
aleSolution = []; 

parameters = [];


if isfile(mbdSolution_filename)
    mdbSolution = csvread(mbdSolution_filename);
    parameters = csvread('simulationParameters.txt');
else
    disp("Cound not generate any visualization! No output file available!");
    return;
end


if isfile(aleSolution_filename)
    aleSolution = csvread(aleSolution_filename);
else
    disp("Cound not generate any visualization! No output file available!");
    return;
end


% Simulation parameters
delta_t = parameters(1);
numIterations = size(aleSolution, 1);

Zdiff = aleSolution(1,60) - mbdSolution(1,60);

% adapt results
mbdSolution(:,60) = mbdSolution(:,60) + Zdiff;
mbdSolution(:,57) = mbdSolution(:,57) + Zdiff;
mbdSolution(:,54) = mbdSolution(:,54) + Zdiff;
mbdSolution(:,51) = mbdSolution(:,51) + Zdiff;

mbdSolution(:,48) = mbdSolution(:,48) + Zdiff;
mbdSolution(:,45) = mbdSolution(:,45) + Zdiff;
mbdSolution(:,42) = mbdSolution(:,42) + Zdiff;
mbdSolution(:,39) = mbdSolution(:,39) + Zdiff;

mbdSolution(:,36) = mbdSolution(:,36) + Zdiff;

diff_all = (mbdSolution - aleSolution) ./ mbdSolution;



