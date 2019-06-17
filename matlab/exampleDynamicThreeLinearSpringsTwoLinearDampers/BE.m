clc;
clear;
close all;
format long e;
% Dynamic Backward Euler problem: CoSim and Mono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_1=1;
k_2=2;
k_3=3;
d_1=1/10;
d_2=1/2;
m_1=1e-1;
m_2=2e-1;
m_3=3e-1;
tend=1;
u_init=1;
du_init=0;
n=23; % refinement 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M=[m_1 0 0;
   0 m_2 0;
   0 0 m_3];
D=[d_1 0 0;
   0 d_2 -d_2;
   0 -d_2 d_2];
K=[(k_1+k_2) -k_2 0;
   -k_2 k_2 0;
   0 0 k_3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
t=zeros(2^n+1,1);
u=zeros(2^n+1,1);
v=zeros(2^n+1,1);
w=zeros(2^n+1,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1/(2^n)
systemMono=BEsolver(h,tend,M,D,K,[u_init; 0; 0],[du_init; 0; 0]);
% Time loop
j=2;
tic


for i = h:h:tend
   [phi, dphi, ddphi]=systemMono.doSolve([0;0;0]); 
   % Get solution, dsolution and ddsolution at rPull level
    t(j)=i;   
    u(j)=phi(1);  
    v(j)=phi(2); 
    w(j)=phi(3); 
    systemMono.incStepCounter();
    j=j+1;
end
toc
j
% Avoid memory leak
delete(systemMono);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write ouput
xy=[t u];
save('solutionBEU.dat', 'xy','-ascii');
xy=[t v];
save('solutionBEV.dat', 'xy','-ascii');
xy=[t w];
save('solutionBEW.dat', 'xy','-ascii');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
figure();
plot(t,u);
hold on;
    plot(t,v);
    hold on;
    plot(t,w);
    legend('u','v','w');