clear; clc; close all;
format long e;
%%
a = 1;
b = 2;
c = 3;

%%
global l_long_fl;
l_long_fl = sqrt(6);
global l_long_fr;
l_long_fr = sqrt(10);
global l_long_rl;
l_long_rl = sqrt(7);
global l_long_rr;
l_long_rr = sqrt(11);
global l_lat_fl;
l_lat_fl = sqrt(3);
global l_lat_fr;
l_lat_fr = sqrt(13);
global l_lat_rl;
l_lat_rl = sqrt(5);
global l_lat_rr;
l_lat_rr = sqrt(11);
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
%%
% evaluation density
dl = 0.1;
dw = 0.1;
% calc min and max length to evaluate (just for now)
l_min = 0.1;
l_max = 5;
w_min = -1;
w_max = 1;
% grid size
l_d = round((l_max - l_min)/dl);
w_d = round((w_max - w_min)/dw);
% grid value allocation
global k_grid;
global X;
global Y;
global Z;
[X,Y,Z] = ndgrid(l_min:dl:l_max,w_min:dw:w_max,w_min:dw:w_max);
k_grid = zeros(size(X));
% k response function
k = @(l,p,t)(c*l*l + b*l + a + p + t);
% fill in grid values
for i = 0:l_d
    for p = 0:w_d
        for t = 0:w_d
            k_grid(i+1,p+1,t+1)= k(l_min+i*dl, w_min+p*dw, w_min+t*dw);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = 1000;
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
M=diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K = get_K(u_n);
D = K *0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_pre=(1/(h*h))*M+(1/h)*D;
B=((2/(h*h))*M+(1/h)*D);
f_n_p_1=[1.1; zeros(10,1)];

%% Time loop
j=2;
tic
for i = h:h:tend
    K = get_K(u_n);
    A=A_pre+K;
    u_n_p_1=A\(B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
   % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;   
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
end
toc

% Some plots
figure();
plot(t,u_sol(:,1)); grid on;
legend;


function K = get_K(x)
    global l_long_fl;
    global l_long_fr;
    global l_long_rl;
    global l_long_rr;
    global l_lat_fl;
    global l_lat_fr;
    global l_lat_rl;
    global l_lat_rr;
    long_fl = l_long_fl -x(4);
    lat_fl = l_lat_fl - x(5);
    long_fr = l_long_fr - x(6);
    lat_fr = l_lat_fr - x(7);
    long_rl = l_long_rl - x(8);
    lat_rl = l_lat_rl - x(9);
    long_rr = l_long_rr - x(10);
    lat_rr = l_lat_rr - x(11);
    
    global k_grid;
    global X;
    global Y;
    global Z;
    
    % for spline interp3(X,Y,Z,k_grid,x(2),long_fl,x(3),'spline')
    k_body_fl=interpn(X,Y,Z,k_grid,long_fl,x(2),x(3));
    k_tyre_fl=interpn(X,Y,Z,k_grid,lat_fl,x(2),x(3));
    k_body_fr=interpn(X,Y,Z,k_grid,long_fr,x(2),x(3));
    k_tyre_fr=interpn(X,Y,Z,k_grid,lat_fr,x(2),x(3));
    k_body_rl=interpn(X,Y,Z,k_grid,long_rl,x(2),x(3));
    k_tyre_rl=interpn(X,Y,Z,k_grid,lat_rl,x(2),x(3));
    k_body_rr=interpn(X,Y,Z,k_grid,long_rr,x(2),x(3));
    k_tyre_rr=interpn(X,Y,Z,k_grid,lat_rr,x(2),x(3));
    
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