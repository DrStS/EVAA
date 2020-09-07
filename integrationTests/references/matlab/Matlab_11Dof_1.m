%%
% Copyright &copy; 2020, Nicola Zimmermann, Munich \n
% All rights reserved. \n
% 
% This file is part of EVAA.
% 
% EVAA is free software: you can redistribute it and/or modify \n
% it under the terms of the GNU General Public  License as published by \n
% the Free Software Foundation, either version 3 of the License, or \n
% (at your option) any later version. \n
% 
% EVAA is distributed in the hope that it will be useful, \n
% but WITHOUT ANY WARRANTY; without even the implied warranty of \n
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
% GNU General Public License for more details. \n
% 
% You should have received a copy of the GNU General Public License \n
% along with EVAA.  If not, see <a
% href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
% 
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Intel Math Kernel Libraries(MKL) (or a modified version of that
% library), containing parts covered by the terms of the license of the MKL,
% the licensors of this Program grant you additional permission to convey the
% resulting work.
% 
% DESCRIPTION
% This file is used to compute displacements/ rotations
% of the 11-Dof 2-track-model
% with a backward Euler time integration scheme 
% for linear spring stiffness

%%
% clc;
clear;
close all;

%%
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
k_stab_f=20e3;
k_stab_r=4e3;
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
no_timesteps=1000;
u_init=0;
du_init=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M=diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K=[k_body_fl+k_body_fr+k_body_rl+k_body_rr, -k_body_fl*l_lat_fl+k_body_fr*l_lat_fr-k_body_rl*l_lat_rl+k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0  l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, +l_long_fl*l_lat_fl*k_body_fl-l_lat_fr*l_long_fr*k_body_fr-l_long_rl*l_lat_rl*k_body_rl+l_long_rr*l_lat_rr*k_body_rr, +l_lat_fl*k_body_fl, 0, -l_lat_fr*k_body_fr, 0, +l_lat_rl*k_body_rl, 0, -l_lat_rr*k_body_rr, 0;
   0 0 l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0 l_long_fr*k_body_fr 0 -l_long_rl*k_body_rl 0 -l_long_rr*k_body_rr 0; 
   0 0 0 k_body_fl+k_tyre_fl -k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 0 k_body_fr+k_tyre_fr -k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 0 k_body_rl+k_tyre_rl -k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 0 k_body_rr+k_tyre_rr -k_tyre_rr;
   0 0 0 0 0 0 0 0 0 0 k_tyre_rr];
K_stab=[0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 +k_stab_f 0 -k_stab_f 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 +k_stab_f 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 +k_stab_r 0 -k_stab_r 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 +k_stab_r 0;
    0 0 0 0 0 0 0 0 0 0 0];
K=K+K_stab;
K=K+K'-diag(diag(K));
D = K *0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=tend/no_timesteps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
dim_system = length(M);
% t=zeros(no_timesteps+1,1);
u_sol=zeros(1000+1,dim_system);
u_n_p_1=zeros(dim_system,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_n=[u_init; zeros(dim_system-1,1)];
u_n_m_1=u_n-h*[du_init; zeros(dim_system-1,1)];
A=((1/(h*h))*M+(1/h)*D+K);
B=((2/(h*h))*M+(1/h)*D);
f_n_p_1=[1.1e3; zeros(dim_system-1,1)];
t = 0:h:(no_timesteps)*h;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%displacement boundary condition
u_f=zeros(11,no_timesteps+1);
u_idx=[5, 7, 9, 11];
% u_idx=[];
% u_f(u_idx,:)=0.1*ones(4,1)*sin(t);
u_f(u_idx,:)=0.01*ones(length(u_idx),1)*t;
% u_f(u_idx,:)=0.0*ones(4,length(t));


%%
% Time loop
% j=2;
tic
for j = 2:length(t)
    rhs = (B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1-A*u_f(:,j));
    A_m=A;
    A_m(u_idx,:)=zeros(length(u_idx),11);
    A_m(:,u_idx)=zeros(11,length(u_idx));
    A_m(u_idx,u_idx)=eye(length(u_idx));
    rhs(u_idx)=zeros(length(u_idx),1);
    rhs=rhs+u_f(:,j);
    u_n_p_1 = A_m\rhs;
%     t(j)=i;   
    u_sol(j,:) = u_n_p_1;
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
%     j = j+1;
end
toc
final_displacement = u_n_p_1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
figure;
subplot(1,5,1);
plot(t,u_sol(:,1)); grid on; legend('z_{CG}'); 
subplot(1,5,2); 
plot(t,u_sol(:,2)); grid on; legend('r_x');
subplot(1,5,3);
plot(t,u_sol(:,3)); grid on; legend('r_y');
subplot(1,5,4);
plot(t,u_sol(:,4)); grid on; legend('z_2');
subplot(1,5,5);
plot(t,u_sol(:,end)); grid on; legend('z_9');

disp(u_sol(end,1:3))

%%
writematrix([t',u_sol],'Matlab_11Dof_1u.dat','Delimiter',',');