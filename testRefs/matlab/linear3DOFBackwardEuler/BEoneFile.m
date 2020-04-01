clc;
clear;
close all;
format long e;
% Dynamic Backward Euler problem
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
u_n_p_1=zeros(3,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1/(2^n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_n=[u_init; 0; 0];
u_n_m_1=[u_init; 0; 0]-h*[du_init; 0; 0];
A=((1/(h*h))*M+(1/h)*D+K);
B=((2/(h*h))*M+(1/h)*D);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop
j=2;
tic
for i = h:h:tend
    u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1);
   % Get solution
    t(j)=i;   
    u(j)=u_n_p_1(1);  
    v(j)=u_n_p_1(2); 
    w(j)=u_n_p_1(3); 
   
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
figure();
plot(t,u);
hold on;
plot(t,v);
hold on;
plot(t,w);
legend('u','v','w');