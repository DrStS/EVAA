clc;
clear all;
close all;
format long e;
% Dynamic Backward Euler problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_body_fl=28e3*0.69;
k_tyre_fl=260e3;%260
k_body_fr=28e3*0.69;
k_tyre_fr=260e3;
k_body_rl=16e3*0.82;
k_tyre_rl=260e3;
k_body_rr=16e3*0.82;
k_tyre_rr=260e3;
l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=2*0.8458;
l_lat_fr=2*0.8458;
l_lat_rl=2*0.84;
l_lat_rr=2*0.84;
mass_Body=1936;
I_body_xx=640;
I_body_yy=4800;
mass_wheel_fl=145/2;
mass_tyre_fl=0;
mass_wheel_fr=145/2;
mass_tyre_fr=0;
mass_wheel_rl=135/2;
mass_tyre_rl=0;
mass_wheel_rr=135/2;
mass_tyre_rr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tend=1;
u_init=0;
du_init=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M=diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K=[k_body_fl+k_body_fr+k_body_rl+k_body_rr, k_body_fl*l_lat_fl-k_body_fr*l_lat_fr+k_body_rl*l_lat_rl-k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0  l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, -l_long_fl*l_lat_fl*k_body_fl+l_lat_fr*l_long_fr*k_body_fr+l_long_rl*l_lat_rl*k_body_rl-l_long_rr*l_lat_rr*k_body_rr, -l_lat_fl*k_body_fl, 0, l_lat_fr*k_body_fr, 0, -l_lat_rl*k_body_rl, 0, l_lat_rr*k_body_rr, 0;
   0 0 l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0 l_long_fr*k_body_fr 0 -l_long_rl*k_body_rl 0 -l_long_rr*k_body_rr 0; 
   0 0 0 k_body_fl+k_tyre_fl -k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 0 k_body_fr+k_tyre_fr -k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 0 k_body_rl+k_tyre_rl -k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 0 k_body_rr+k_tyre_rr -k_tyre_rr;
   0 0 0 0 0 0 0 0 0 0 k_tyre_rr];
K=K+K'-diag(diag(K));
D = K *0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced System due to Dirichlet BC (fixed to road)
Kred = K;
Kred(:,11) = [];
Kred(11,:) = [];
Kred(:,9) = [];
Kred(9,:) = [];
Kred(:,7) = [];
Kred(7,:) = [];
Kred(:,5) = [];
Kred(5,:) = [];
Dred = D;
Dred(:,11) = [];
Dred(11,:) = [];
Dred(:,9) = [];
Dred(9,:) = [];
Dred(:,7) = [];
Dred(7,:) = [];
Dred(:,5) = [];
Dred(5,:) = [];
Mred = M;
Mred(:,11) = [];
Mred(11,:) = [];
Mred(:,9) = [];
Mred(9,:) = [];
Mred(:,7) = [];
Mred(7,:) = [];
Mred(:,5) = [];
Mred(5,:) = [];

% K(:,11) = [];
% K(11,:) = [];
% K(:,9) = [];
% K(9,:) = [];
% K(:,7) = [];
% K(7,:) = [];
% K(:,5) = [];
% K(5,:) = [];
% D(:,11) = [];
% D(11,:) = [];
% D(:,9) = [];
% D(9,:) = [];
% D(:,7) = [];
% D(7,:) = [];
% D(:,5) = [];
% D(5,:) = [];
% M(:,11) = [];
% M(11,:) = [];
% M(:,9) = [];
% M(9,:) = [];
% M(:,7) = [];
% M(7,:) = [];
% M(:,5) = [];
% M(5,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
dim_system = length(M);
t=zeros(1000+1,1);
u_sol=zeros(1000+1,dim_system);
u_n_p_1=zeros(dim_system,1);
u_sol_red=zeros(1000+1,7);
u_n_p_1_red=zeros(7,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1/(1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_n=[u_init; zeros(dim_system-1,1)];
u_n_m_1=[u_init; zeros(dim_system-1,1)]-h*[du_init; zeros(dim_system-1,1)];
u_n_red=[u_init; zeros(6,1)];
u_n_m_1_red=[u_init; zeros(6,1)]-h*[du_init; zeros(6,1)];
A=((1/(h*h))*M+(1/h)*D+K);
B=((2/(h*h))*M+(1/h)*D);
Ared=((1/(h*h))*Mred+(1/h)*Dred+Kred);
Bred=((2/(h*h))*Mred+(1/h)*Dred);
f_n_p_1=[1.1e3; zeros(dim_system-1,1)];
f_n_p_1_red=[1.1e3; zeros(6,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static solution
u_stat=Kred\f_n_p_1(1:7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigenfrequency solution with bc
eigFreq=sqrt(eig(Kred,Mred))/2/pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop
j=2;
idx = [1,2,3,4,6,8,10];
tic
for i = h:h:tend
    u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
    u_n_p_1_red=Ared\(Bred*u_n_red-((1/(h*h))*Mred)*u_n_m_1_red+f_n_p_1_red);
   % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    u_sol_red(j,:)=u_n_p_1_red;
    u_n_m_1_red=u_n_red;
    u_n_red    =u_n_p_1_red;
    j=j+1;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
figure();
plot(t,u_sol(:,1)); grid on;
%plot(t,u_sol_red(:,1)); grid on;
legend;
disp(u_sol(1001,1:11))
disp(u_sol_red(1001,1:3))