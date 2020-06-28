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
% for nonlinear spring stiffness, 
% with parameters set to constant stiffness
% compared to BEoneFile03:
%  - nonlinear stiffness
%%
clear; 
%clc;
%close all;
format long e;
%% stiffness params (quadratic)
a = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3] ;
% b = 1.59;
% c = 1.2;
b=0;
c=0;

%%
global l_long_fl;
global l_long_fr;
global l_long_rl;
global l_long_rr;
global l_lat_fl;
global l_lat_fr;
global l_lat_rl;
global l_lat_rr;
l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=1.6916;
l_lat_fr=1.6916;
l_lat_rl=1.68;
l_lat_rr=1.68;
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
interpolation = "linear";
tend=1;
u_init=0;
du_init=0;
global Corners;
Corners = [ l_long_fl, l_long_fr, -l_long_rl, -l_long_rr;
    l_lat_fl, -l_lat_fr, l_lat_rl, -l_lat_rr;
    0, 0, 0, 0];
%%
% calc min and max length to evaluate (just for now)
l_min = 0.05;
l_max = 0.8;
% grid size
global size;
size = 1000;
% grid value allocation
global k_grid;
global X;
X = zeros(size, 1);
k_grid = zeros(8*size,1);
dl = (l_max-l_min)/(size-1);
% k response function
k = @(l,a)(c*l*l + b*l + a);
% fill in grid values
for i = 0:size-1
    X(i+1) = l_min+i*dl;
end
for i = 1:8
    for j = 1:size
        k_grid((i-1)*size + j)= k(X(j), a(i));
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = 1000;
% h=1/(steps-1);
h=1/steps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
t=zeros(steps+1,1);
u_sol=zeros(steps+1,11);
u_n_p_1=zeros(11,1);
u_n=[u_init; zeros(10,1)];
u_n_m_1=[u_init; zeros(10,1)]-h*[du_init; zeros(10,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K = get_K(u_n);
D = K *0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_pre=(1/(h*h))*M+(1/h)*D;
B=((2/(h*h))*M+(1/h)*D);
f_n_p_1=[1.1e3; zeros(10,1)];


%% Time loop
j=2;
tic
for i = h:h:tend
    K = get_K(u_n);
    A=A_pre+K;
    rhs = B*u_n-((1/(h*h))*M)*u_n_m_1 + f_n_p_1;
    u_n_p_1=A\rhs;
   % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;   
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j = j + 1;
end
toc

% Some plots
% figure();
% plot(t,u_sol(:,1)); grid on;
% legend;
figure;
disp(u_sol(1001,1:3))
subplot(1,3,1);
plot(t,u_sol(:,1)); grid on; legend('u_z'); 
subplot(1,3,2); 
plot(t,u_sol(:,2)); grid on; legend('r_x');
subplot(1,3,3);
plot(t,u_sol(:,3)); grid on; legend('r_y');



function K = get_K(x)
    global l_long_fl;
    global l_long_fr;
    global l_long_rl;
    global l_long_rr;
    global l_lat_fl;
    global l_lat_fr;
    global l_lat_rl;
    global l_lat_rr;
    car = x(1);
    upper_fl = x(4);
    lower_fl = x(5);
    upper_fr = x(6);
    lower_fr = x(7);
    upper_rl = x(8);
    lower_rl = x(9);
    upper_rr = x(10);
    lower_rr = x(11);
    
    global k_grid;
    global X;
    
    global Corners size;
    R = get_R(x);
    currentCorners = R*Corners;
    
    % for spline interp3(X,Y,Z,k_grid,x(2),long_fl,x(3),'spline')
    k_body_fl=interp1(X,k_grid(1:size),0.5 + currentCorners(3,1) - upper_fl + car);
    k_tyre_fl=interp1(X,k_grid(size+1:2*size),0.5 + upper_fl - lower_fl);
    k_body_fr=interp1(X,k_grid(2*size+1:3*size),0.5 + currentCorners(3,2) - upper_fr + car);
    k_tyre_fr=interp1(X,k_grid(3*size+1:4*size),0.5 + upper_fr - lower_fr);
    k_body_rl=interp1(X,k_grid(4*size+1:5*size),0.5 + currentCorners(3,3) - upper_rl + car);
    k_tyre_rl=interp1(X,k_grid(5*size+1:6*size),0.5 + upper_rl -lower_rl);
    k_body_rr=interp1(X,k_grid(6*size+1:7*size),0.5 + currentCorners(3,4) - upper_rr + car);
    k_tyre_rr=interp1(X,k_grid(7*size+1:8*size),0.5 + upper_rr - lower_rr);
    
    K = [k_body_fl+k_body_fr+k_body_rl+k_body_rr, k_body_fl*l_lat_fl-k_body_fr*l_lat_fr+k_body_rl*l_lat_rl-k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0, l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, -l_long_fl*l_lat_fl*k_body_fl+l_lat_fr*l_long_fr*k_body_fr+l_long_rl*l_lat_rl*k_body_rl-l_long_rr*l_lat_rr*k_body_rr, -l_lat_fl*k_body_fl, 0, l_lat_fr*k_body_fr, 0, -l_lat_rl*k_body_rl, 0, l_lat_rr*k_body_rr, 0;
   0, 0, l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0, l_long_fr*k_body_fr, 0, -l_long_rl*k_body_rl, 0, -l_long_rr*k_body_rr, 0; 
   0, 0, 0, k_body_fl+k_tyre_fl, -k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k_body_fr+k_tyre_fr, -k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k_body_rl+k_tyre_rl, -k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k_body_rr+k_tyre_rr, -k_tyre_rr;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rr];
    K=K+K'-diag(diag(K));
end
function R = get_R (x)
 xx = x(2);
 yy = x(3);
 R = [cos(yy), sin(yy)*sin(xx), sin(yy)*cos(xx);
     0, cos(xx), -sin(xx);
     -sin(yy), cos(yy)*sin(xx), cos(yy)*cos(xx)
     ];
end