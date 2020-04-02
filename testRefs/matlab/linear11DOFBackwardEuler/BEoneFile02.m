%%
clc;
clear all;
close all;
format long e;
% Dynamic Backward Euler problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_body_fl=sqrt(2);
k_tyre_fl=sqrt(40);
k_body_fr=sqrt(3);
k_tyre_fr=sqrt(50);
k_body_rl=sqrt(5);
k_tyre_rl=sqrt(60);
k_body_rr=sqrt(7);
k_tyre_rr=sqrt(80);
l_long_fl=sqrt(6);
l_long_fr=sqrt(10);
l_long_rl=sqrt(7);
l_long_rr=sqrt(11);
l_lat_fl=sqrt(3);
l_lat_fr=sqrt(13);
l_lat_rl=sqrt(5);
l_lat_rr=sqrt(11);
mass_Body=sqrt(3e-1);
I_body_xx=sqrt(5e-1);
I_body_yy=sqrt(6e-1);
mass_wheel_fl=sqrt(5e-1);
mass_tyre_fl=sqrt(7e-1);
mass_wheel_fr=sqrt(9e-1);
mass_tyre_fr=sqrt(11e-1);
mass_wheel_rl=sqrt(13e-1);
mass_tyre_rl=sqrt(15e-1);
mass_wheel_rr=sqrt(17e-1);
mass_tyre_rr=sqrt(23e-1);
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
% Allocate memory
t=zeros(1000+1,1);
u_sol=zeros(1000+1,11);
u_n_p_1=zeros(11,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=1/(1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_n=[u_init; zeros(10,1)];
u_n_m_1=[u_init; zeros(10,1)]-h*[du_init; zeros(10,1)];
A=((1/(h*h))*M+(1/h)*D+K);
B=((2/(h*h))*M+(1/h)*D);
f_n_p_1=[1.1; zeros(10,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time loop
j=2;
tic
for i = h:h:tend
    u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
   % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;   
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
end
toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
figure();
plot(t,u_sol(:,1)); grid on;
legend;

%disp(u_sol(501,1))
disp(u_sol(1001,1:3))

