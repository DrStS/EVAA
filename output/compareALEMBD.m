%
% Only compare two results if the initial spring lengths is the same
%


mbdSolution_filename = "mbdSolution.txt";
mbdSolution = [];

aleSolution_filename = "aleSolution.txt";
aleSolution = []; 

parameters = [];


if isfile(mbdSolution_filename)
    mbdSolution = csvread(mbdSolution_filename);
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

Zdiff = aleSolution(2,60) - mbdSolution(2,60);

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

diff_all = (mbdSolution - aleSolution);

index = 35;

normed_pos = vecnorm(diff_all(2:end,index:end),2,2);

t = delta_t * (1:size(diff_all(2:end,:), 1)); 

%compute the DFT
mbdDFT = fft(mbdSolution(2:end,36) - mean(mbdSolution(2:end,36)));
aleDFT = fft(aleSolution(2:end,36) - mean(aleSolution(2:end,36)));
sampleFrequency = 1 / delta_t;
L = length(aleSolution(2:end,36));

P2 = abs(mbdDFT/L);
mbdDFT = P2(1:L/2+1);
mbdDFT(2:end-1) = 2*mbdDFT(2:end-1);


P2 = abs(aleDFT/L);
aleDFT = P2(1:L/2+1);
aleDFT(2:end-1) = 2*aleDFT(2:end-1);

frequency = sampleFrequency * (0:(L/2))/L;


close all
figure
plot(t,diff_all(2:end,index + 0))
hold on
plot(t,diff_all(2:end,index + 2))
plot(t,diff_all(2:end,index + 1))
legend('X', 'Y', 'Z')
title('Difference of the position of the center of gravity')
ylabel('$x_{MBD} - x_{ALE}$','Interpreter','latex','FontSize',16)
xlabel('time [s]');

figure
plot(t,normed_pos)
hold on
legend('Implicit Euler')
title('Norm of the difference of all positions(CoG, wheels tyres)')
ylabel('$||x_{MBD} - x_{ALE}||_2$','Interpreter','latex','FontSize',16)
xlabel('time [s]');

figure
plot(t,mbdSolution(2:end,36),'g')
hold on
plot(t,aleSolution(2:end,36),'r')
legend('MBD','ALE')
title('Evolution of the Z-component of the center of gravity')
ylabel('Z-position [m]')
xlabel('time [s]');

figure
plot(t,diff_all(2:end,index + 0))
hold on
plot(t,diff_all(2:end,index + 2))
plot(t,mbdSolution(2:end,36)-mbdSolution(1,36))
legend('X','Y','Z')
title('Deviation from the trajectory for MBD, Ixx=6400, Iyy=Izz=6400')
ylabel('difference [m]')
xlabel('time [s]')
axis([0 50 -0.08 0.01])

figure
plot(t,zeros(length(t),1))  % zero due to the Lagrangiabn formulation
hold on
plot(t,zeros(length(t),1))  % zero due to the Lagrangiabn formulation
plot(t,aleSolution(2:end,36)-aleSolution(2,36))
legend('X','Y','Z')
title('Deviation from the trajectory for ALE, Ixx=6400, Iyy=Izz=4800')
ylabel('difference [m]')
xlabel('time [s]')
axis([0 50 -0.08 0.01])

figure
plot(frequency(1:end), aleDFT(1:end))
hold on
plot(frequency(1:end), mbdDFT(1:end))
legend('MBD','ALE')
title('Amplitude spectrum of the Z-component of the center of gravity')
xlabel('f [Hz]')
ylabel('|P1(f)|')






