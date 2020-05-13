clear; clc; close all;
%% make lookup
ak = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3];
ad = ak/10;
b = 1.59;
c = 1.2;
%b = 0;
%c = 0;
%% zero lengths

L1 = 0.2;
L2 = 0.2;
L3 = 0.2;
L4 = 0.2;
L5 = 0.2;
L6 = 0.2;
L7 = 0.2;
L8 = 0.2;
%%
global k_grid d_grid size_grid l_min l_max dl;
% calc min and max length to evaluate (just for now)
l_min = 0.01;
l_max = 2;
% grid size
size_grid = 1001;
% grid value allocation
global X;
X = zeros(size_grid, 1);
k_grid = zeros(8*size_grid,1);
d_grid = zeros(8*size_grid,1);
dl = (l_max-l_min)/(size_grid-1);
% k response function
k = @(l,a)(c*l*l + b*l + a);
%k = @(l,a)(a);
% fill in grid values
for i = 0:size_grid-1
    X(i+1) = l_min+i*dl;
end
for i = 1:8 
    for j = 1:size_grid
        k_grid((i-1)*size_grid + j)= k(X(j), ak(i));
    end
end
for i = 1:8 
    for j = 1:size_grid
        d_grid((i-1)*size_grid + j)= k(X(j), ad(i));
    end
end
%dlmwrite('grid.txt',[X,k_grid(1:size_grid)],'delimiter',',','precision',9);

global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 
global d_spline1 d_spline2 d_spline3 d_spline4 d_spline5 d_spline6 d_spline7 d_spline8;
global d_der1 d_der2 d_der3 d_der4 d_der5 d_der6 d_der7 d_der8; 
k_spline1 = spline(X,k_grid(1:size_grid));
k_der1 = fnder(k_spline1,1);
k_spline2 = spline(X,k_grid(size_grid+1:2*size_grid));
k_der2 = fnder(k_spline2,1);
k_spline3 = spline(X,k_grid(2*size_grid+1:3*size_grid));
k_der3 = fnder(k_spline3,1);
k_spline4 = spline(X,k_grid(3*size_grid+1:4*size_grid));
k_der4 = fnder(k_spline4,1);
k_spline5 = spline(X,k_grid(4*size_grid+1:5*size_grid));
k_der5 = fnder(k_spline5,1);
k_spline6 = spline(X,k_grid(5*size_grid+1:6*size_grid));
k_der6 = fnder(k_spline6,1);
k_spline7 = spline(X,k_grid(6*size_grid+1:7*size_grid));
k_der7 = fnder(k_spline7,1);
k_spline8 = spline(X,k_grid(7*size_grid+1:8*size_grid));
k_der8 = fnder(k_spline8,1);
d_spline1 = spline(X,d_grid(1:size_grid));
d_der1 = fnder(d_spline1,1);
d_spline2 = spline(X,d_grid(size_grid+1:2*size_grid));
d_der2 = fnder(k_spline2,1);
d_spline3 = spline(X,d_grid(2*size_grid+1:3*size_grid));
d_der3 = fnder(d_spline3,1);
d_spline4 = spline(X,d_grid(3*size_grid+1:4*size_grid));
d_der4 = fnder(d_spline4,1);
d_spline5 = spline(X,d_grid(4*size_grid+1:5*size_grid));
d_der5 = fnder(d_spline5,1);
d_spline6 = spline(X,d_grid(5*size_grid+1:6*size_grid));
d_der6 = fnder(d_spline6,1);
d_spline7 = spline(X,d_grid(6*size_grid+1:7*size_grid));
d_der7 = fnder(d_spline7,1);
d_spline8 = spline(X,d_grid(7*size_grid+1:8*size_grid));
d_der8 = fnder(d_spline8,1);

%%
global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr;
global x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11;
global d1 d2 d3 d4 d5 d6 d7 d8
global l1 l2 l3 l4 l5 l6 l7 l8
%% symbolic derivative
syms l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr
k = sym('k',[8 1]); 
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
syms d1 d2 d3 d4 d5 d6 d7 d8
syms l1_sym l2_sym l3_sym l4_sym l5_sym l6_sym l7_sym l8_sym
syms dKdx_x_symb
l1_sym = L1 + x1 + l_lat_fl*x2 - l_long_fl*x3  - x4;
l2_sym = L2 + x4 - x5;
l3_sym = L3 + x1 - l_lat_fr*x2 - l_long_fr*x3 - x6;
l4_sym = L4 + x6 - x7;
l5_sym = L5 + x1 + l_lat_rl*x2 + l_long_rl*x3 - x8;
l6_sym = L6 + x8 - x9;
l7_sym = L7 + x1 - x10 - l_lat_rr*x2 + l_long_rr*x3;
l8_sym = L8 + x10 - x11;

dKdx_x_symb = [(d1 + d3 + d5 + d7)*x1 - d7*x10 + (d1*l_lat_fl - d3*l_lat_fr + ...
d5*l_lat_rl - d7*l_lat_rr)*x2 - d1*x4 - d3*x6 - d5*x8 + x3*(-(d3*l_long_fr)...
+ d5*l_long_rl + d7*l_long_rr - d1*(l_long_fl))	-(d3*l_lat_fr*x1) + ...
d5*l_lat_rl*x1 - d7*l_lat_rr*x1 + d7*l_lat_rr*x10 + d3*l_lat_fr^2*x2 + ...
d5*l_lat_rl^2*x2 + d7*l_lat_rr^2*x2 + d3*l_lat_fr*l_long_fr*x3 + ...
d5*l_lat_rl*l_long_rl*x3 - d7*l_lat_rr*l_long_rr*x3 + d1*l_lat_fl*(x1 + ...
l_lat_fl*x2 - x4) + d3*l_lat_fr*x6 - d5*l_lat_rl*x8 - ...
d1*l_lat_fl*x3*(l_long_fl)	d5*l_long_rl*x1 + d7*l_long_rr*x1 - ...
d7*l_long_rr*x10 + d5*l_lat_rl*l_long_rl*x2 - d7*l_lat_rr*l_long_rr*x2 + ...
d5*l_long_rl^2*x3 + d7*l_long_rr^2*x3 + d3*l_long_fr*(-x1 + l_lat_fr*x2 + ...
l_long_fr*x3 + x6) - d5*l_long_rl*x8 - d1*(x1 + l_lat_fl*x2 - ...
x4)*(l_long_fl) + d1*x3*(l_long_fl)^2	d1*(-x1 - l_lat_fl*x2 + x4 + ...
x3*(l_long_fl))	0	d3*(-x1 + l_lat_fr*x2 + l_long_fr*x3 + x6)	0	-(d5*(x1 + ...
l_lat_rl*x2 + l_long_rl*x3 - x8))	0	d7*(-x1 + x10 + l_lat_rr*x2 - ...
l_long_rr*x3)	0;-(d3*l_lat_fr*x1) + d5*l_lat_rl*x1 - d7*l_lat_rr*x1 + ...
d7*l_lat_rr*x10 + d3*l_lat_fr^2*x2 + d5*l_lat_rl^2*x2 + d7*l_lat_rr^2*x2 + ...
d3*l_lat_fr*l_long_fr*x3 + d5*l_lat_rl*l_long_rl*x3 - d7*l_lat_rr*l_long_rr*x3 + ...
d1*l_lat_fl*(x1 + l_lat_fl*x2 - x4) + d3*l_lat_fr*x6 - d5*l_lat_rl*x8 - ...
d1*l_lat_fl*x3*(l_long_fl)	d3*l_lat_fr^2*x1 + d5*l_lat_rl^2*x1 + ...
d7*l_lat_rr^2*x1 - d7*l_lat_rr^2*x10 - d3*l_lat_fr^3*x2 + d5*l_lat_rl^3*x2 - ...
d7*l_lat_rr^3*x2 - d3*l_lat_fr^2*l_long_fr*x3 + d5*l_lat_rl^2*l_long_rl*x3 + ...
d7*l_lat_rr^2*l_long_rr*x3 + d1*l_lat_fl^2*(x1 + l_lat_fl*x2 - x4) - ...
d3*l_lat_fr^2*x6 - d5*l_lat_rl^2*x8 - d1*l_lat_fl^2*x3*(l_long_fl)	...
d3*l_lat_fr*l_long_fr*x1 + d5*l_lat_rl*l_long_rl*x1 - d7*l_lat_rr*l_long_rr*x1 + ...
d7*l_lat_rr*l_long_rr*x10 - d3*l_lat_fr^2*l_long_fr*x2 + ...
d5*l_lat_rl^2*l_long_rl*x2 + d7*l_lat_rr^2*l_long_rr*x2 - ...
d3*l_lat_fr*l_long_fr^2*x3 + d5*l_lat_rl*l_long_rl^2*x3 - ...
d7*l_lat_rr*l_long_rr^2*x3 - d1*l_lat_fl*l_long_fl*(x1 + l_lat_fl*x2 - x4) - ...
d3*l_lat_fr*l_long_fr*x6 - d5*l_lat_rl*l_long_rl*x8 + ...
d1*l_lat_fl*l_long_fl*x3*(l_long_fl)	d1*l_lat_fl*(-x1 - l_lat_fl*x2 + x4 + ...
x3*(l_long_fl))	0	-(d3*l_lat_fr*(-x1 + l_lat_fr*x2 + l_long_fr*x3 + x6))	0	...
-(d5*l_lat_rl*(x1 + l_lat_rl*x2 + l_long_rl*x3 - x8))	0	d7*l_lat_rr*(x1 - x10 ...
- l_lat_rr*x2 + l_long_rr*x3)	0;d5*l_long_rl*x1 + d7*l_long_rr*x1 - ...
d7*l_long_rr*x10 + d5*l_lat_rl*l_long_rl*x2 - d7*l_lat_rr*l_long_rr*x2 + ...
d5*l_long_rl^2*x3 + d7*l_long_rr^2*x3 + d3*l_long_fr*(-x1 + l_lat_fr*x2 + ...
l_long_fr*x3 + x6) - d5*l_long_rl*x8 - d1*(x1 + l_lat_fl*x2 - ...
x4)*(l_long_fl) + d1*x3*(l_long_fl)^2	d3*l_lat_fr*l_long_fr*x1 + ...
d5*l_lat_rl*l_long_rl*x1 - d7*l_lat_rr*l_long_rr*x1 + d7*l_lat_rr*l_long_rr*x10 - ...
d3*l_lat_fr^2*l_long_fr*x2 + d5*l_lat_rl^2*l_long_rl*x2 + ...
d7*l_lat_rr^2*l_long_rr*x2 - d3*l_lat_fr*l_long_fr^2*x3 + ...
d5*l_lat_rl*l_long_rl^2*x3 - d7*l_lat_rr*l_long_rr^2*x3 - ...
d1*l_lat_fl*l_long_fl*(x1 + l_lat_fl*x2 - x4) - d3*l_lat_fr*l_long_fr*x6 - ...
d5*l_lat_rl*l_long_rl*x8 + d1*l_lat_fl*l_long_fl*x3*(l_long_fl)	...
d3*l_long_fr^2*x1 + d5*l_long_rl^2*x1 + d7*l_long_rr^2*x1 - ...
d7*l_long_rr^2*x10 - d3*l_lat_fr*l_long_fr^2*x2 + d5*l_lat_rl*l_long_rl^2*x2 - ...
d7*l_lat_rr*l_long_rr^2*x2 - d3*l_long_fr^3*x3 + d5*l_long_rl^3*x3 + ...
d7*l_long_rr^3*x3 + d1*l_long_fl^2*(x1 + l_lat_fl*x2 - x4) - ...
d3*l_long_fr^2*x6 - d5*l_long_rl^2*x8 - d1*l_long_fl^2*x3*(l_long_fl)	...
d1*l_long_fl*(x1 + l_lat_fl*x2 - x4 - x3*(l_long_fl))	0	-(d3*l_long_fr*(-x1 ...
+ l_lat_fr*x2 + l_long_fr*x3 + x6))	0	-(d5*l_long_rl*(x1 + l_lat_rl*x2 + ...
l_long_rl*x3 - x8))	0	d7*l_long_rr*(-x1 + x10 + l_lat_rr*x2 - l_long_rr*x3)	...
0;d1*(-x1 - l_lat_fl*x2 + x4 + x3*(l_long_fl))	d1*l_lat_fl*(-x1 - ...
l_lat_fl*x2 + x4 + x3*(l_long_fl))	d1*l_long_fl*(x1 + l_lat_fl*x2 - x4 - ...
x3*(l_long_fl))	d1*(x1 + l_lat_fl*x2 - x4) + d2*(x4 - x5) - ...
d1*x3*(l_long_fl)	d2*(-x4 + x5)	0	0	0	0	0	0;0	0	0	d2*(-x4 + x5)	d2*(x4 ...
- x5)	0	0	0	0	0	0;d3*(-x1 + l_lat_fr*x2 + l_long_fr*x3 + x6)	...
-(d3*l_lat_fr*(-x1 + l_lat_fr*x2 + l_long_fr*x3 + x6))	-(d3*l_long_fr*(-x1 + ...
l_lat_fr*x2 + l_long_fr*x3 + x6))	0	0	d3*(x1 - l_lat_fr*x2 - l_long_fr*x3 - ...
x6) + d4*(x6 - x7)	d4*(-x6 + x7)	0	0	0	0;0	0	0	0	0	d4*(-x6 + x7)	...
d4*(x6 - x7)	0	0	0	0;-(d5*(x1 + l_lat_rl*x2 + l_long_rl*x3 - x8))	...
-(d5*l_lat_rl*(x1 + l_lat_rl*x2 + l_long_rl*x3 - x8))	-(d5*l_long_rl*(x1 + ...
l_lat_rl*x2 + l_long_rl*x3 - x8))	0	0	0	0	d5*(x1 + l_lat_rl*x2 + l_long_rl*x3 ...
- x8) + d6*(x8 - x9)	d6*(-x8 + x9)	0	0;0	0	0	0	0	0	0	d6*(-x8 + x9)	...
d6*(x8 - x9)	0	0;d7*(-x1 + x10 + l_lat_rr*x2 - l_long_rr*x3)	...
d7*l_lat_rr*(x1 - x10 - l_lat_rr*x2 + l_long_rr*x3)	d7*l_long_rr*(-x1 + x10 + ...
l_lat_rr*x2 - l_long_rr*x3)	0	0	0	0	0	0	d8*(x10 - x11) + d7*(x1 - x10 - ...
l_lat_rr*x2 + l_long_rr*x3)	d8*(-x10 + x11);0	0	0	0	0	0	0	0	0	d8*(-x10 + ...
x11)	d8*(x10 - x11)];


%% parameters
% time
num_iter = 1e3;
delta_t = 1e-3; 
t = 0:delta_t:(num_iter - 1)*delta_t;

tol = 1e-7;
y = zeros(length(t),11);
order = zeros(length(t),1);
err_arr = zeros(length(t),1);

l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=1.6916;
l_lat_fr=1.6916;
l_lat_rl=1.68;
l_lat_rr=1.68;
mass_Body=1936;
I_body_xx=6400;
I_body_yy=4800;
mass_wheel_fl=72.5;
mass_tyre_fl=30;
mass_wheel_fr=72.5;
mass_tyre_fr=30;
mass_wheel_rl=67.5;
mass_tyre_rl=30;
mass_wheel_rr=67.5;
mass_tyre_rr=30;
u_n = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
%% init values
x1 = u_n(1);
x2 = u_n(2);
x3 = u_n(3);
x4 = u_n(4);
x5 = u_n(5);
x6 = u_n(6);
x7 = u_n(7);
x8 = u_n(8);
x9 = u_n(9);
x10 = u_n(10);
x11 = u_n(11);
l1 = eval(l1_sym);
l2 = eval(l2_sym);
l3 = eval(l3_sym);
l4 = eval(l4_sym);
l5 = eval(l5_sym);
l6 = eval(l6_sym);
l7 = eval(l7_sym);
l8 = eval(l8_sym);

u_n_m_1 = u_n;
u_n_m_1(1)=u_n(1);
u_n_p_1=u_n;

M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
rhs =[1.1e3; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
%rhs = diag(M)*(-9.81);
rhs(2:3) = 0;
M_div_h2 = M / (delta_t * delta_t);
%%
fixed = true;
euler = true;
bdf2 = true;
interpolation = false;
if bdf2
    B = (6)*M_div_h2;
    C = (-11/2)*M_div_h2;
    D = (2)*M_div_h2;
    E = (-1/4)*M_div_h2;
    u_n_m_3 = zeros(11,1);
    u_n_m_2 = zeros(11,1);
end
%% time steps
idx = [5,7,9,11];
f_newton_BE = @(y_curr,y1,y2,K,rhs)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);
f_newton_Bdf2 = @(y_curr,y1,y2,y3,y4,K,rhs)( ((9/4)*M_div_h2 + K) * y_curr - B * y1 - C * y2 - D * y3 - E * y4 - rhs);
%dKcols_dk = eval(dKcols_dk); % evaluate it numerically
condition = [];
d = 0;
%write to output csv
iterationsRES = zeros(length(t),2);
if (~interpolation)
    k_body_fl=28e3*0.69;
    k_tyre_fl=260e3;%260
    k_body_fr=28e3*0.69;
    k_tyre_fr=260e3;
    k_body_rl=16e3*0.82;
    k_tyre_rl=260e3;
    k_body_rr=16e3*0.82;
    k_tyre_rr=260e3;
    
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
end
%K = get_K();
tic
for i = 1: length(t)
    if interpolation
        K = get_K();
    end
    % J += -(K(tidx,tidx)+ dKdx(tidx,tidx)*u_n_p_1)
    if fixed
        newForce = K * u_n_p_1;
        rhs(5) = newForce(5);
        rhs(7) = newForce(7);
        rhs(9) = newForce(9);
        rhs(11) = newForce(11);
        if euler
            u_n_p_1 = ( M_div_h2 + K )\ (2 * M_div_h2 * u_n - M_div_h2 * u_n_m_1 + rhs);
        elseif bdf2
            u_n_p_1 = ((9/4)*M_div_h2 + K)\ (B*u_n + C*u_n_m_1 + D*u_n_m_2 + E*u_n_m_3 + rhs);
        end
        if interpolation
            x1 = u_n_p_1(1);
            x2 = u_n_p_1(2);
            x3 = u_n_p_1(3);
            x4 = u_n_p_1(4);
            x5 = u_n_p_1(5);
            x6 = u_n_p_1(6);
            x7 = u_n_p_1(7);
            x8 = u_n_p_1(8);
            x9 = u_n_p_1(9);
            x10 = u_n_p_1(10);
            x11 = u_n_p_1(11);
            l1 = eval(l1_sym);
            l2 = eval(l2_sym);
            l3 = eval(l3_sym);
            l4 = eval(l4_sym);
            l5 = eval(l5_sym);
            l6 = eval(l6_sym);
            l7 = eval(l7_sym);
            l8 = eval(l8_sym);
            len = [l1;l2;l3;l4;l5;l6;l7;l8];
            K = get_K();
        end
        newForce = K * u_n_p_1;
        rhs(5) = newForce(5);
        rhs(7) = newForce(7);
        rhs(9) = newForce(9);
        rhs(11) = newForce(11);
    else
        if euler
           u_n_p_1 = ( M_div_h2 + K )\ (2 * M_div_h2 * u_n - M_div_h2 * u_n_m_1 + rhs);
        elseif bdf2
            vec_rhs = B*u_n + C*u_n_m_1 + D*u_n_m_2 + E*u_n_m_3;
            u_n_p_1 = ((9/4)*M_div_h2 + K)\ (vec_rhs + rhs);
        end
        if interpolation
            x1 = u_n_p_1(1);
            x2 = u_n_p_1(2);
            x3 = u_n_p_1(3);
            x4 = u_n_p_1(4);
            x5 = u_n_p_1(5);
            x6 = u_n_p_1(6);
            x7 = u_n_p_1(7);
            x8 = u_n_p_1(8);
            x9 = u_n_p_1(9);
            x10 = u_n_p_1(10);
            x11 = u_n_p_1(11);
            l1 = eval(l1_sym);
            l2 = eval(l2_sym);
            l3 = eval(l3_sym);
            l4 = eval(l4_sym);
            l5 = eval(l5_sym);
            l6 = eval(l6_sym);
            l7 = eval(l7_sym);
            l8 = eval(l8_sym);
            len = [l1;l2;l3;l4;l5;l6;l7;l8];
            K = get_K();
        end
    end
    
    if euler
        r = f_newton_BE(u_n_p_1, u_n, u_n_m_1, K, rhs);
    elseif bdf2
        r = f_newton_Bdf2(u_n_p_1, u_n, u_n_m_1, u_n_m_2, u_n_m_3, K, rhs);
    end
   % i
    iter = 0;
    err = [];
    delta = [];
    init_err = r;
    while 1
        iter = iter + 1;
        temp = [];
        if interpolation
            getdk_dx();
            dKdx_x = eval(dKdx_x_symb);
        end
        if euler
            J = M_div_h2 + K;
            if interpolation
                J = J + dKdx_x;
            end
            if fixed
               J(idx,:) = M_div_h2(idx,:);
            end
        elseif bdf2
            J = (9.0/4.0) * M_div_h2 + K;
            if interpolation
                J = J + dKdx_x;
            end
            if fixed
               J(idx,:) = (9.0/4.0) * M_div_h2(idx,:);
            end
        end
        J
        Delta = -J\r;
        u_n_p_1 = Delta + u_n_p_1; 
        if interpolation
            % update values
            x1 = u_n_p_1(1);
            x2 = u_n_p_1(2);
            x3 = u_n_p_1(3);
            x4 = u_n_p_1(4);
            x5 = u_n_p_1(5);
            x6 = u_n_p_1(6);
            x7 = u_n_p_1(7);
            x8 = u_n_p_1(8);
            x9 = u_n_p_1(9);
            x10 = u_n_p_1(10);
            x11 = u_n_p_1(11);
            l1 = eval(l1_sym);
            l2 = eval(l2_sym);
            l3 = eval(l3_sym);
            l4 = eval(l4_sym);
            l5 = eval(l5_sym);
            l6 = eval(l6_sym);
            l7 = eval(l7_sym);
            l8 = eval(l8_sym);
            len = [l1;l2;l3;l4;l5;l6;l7;l8];
%        for l=1:8
%            if len(l)<l_min
%                disp(len(l))
%                return
%            elseif len(l)>l_max
%                disp(len(l))
%                return
%            end
%        end
            
            K = get_K();
        end
        if fixed
            newForce = K * u_n_p_1;
            rhs(5) = newForce(5);
            rhs(7) = newForce(7);
            rhs(9) = newForce(9);
            rhs(11) = newForce(11); 
        end
        if euler
            r = f_newton_BE(u_n_p_1, u_n, u_n_m_1, K,rhs);
        elseif bdf2
            r = f_newton_Bdf2(u_n_p_1, u_n, u_n_m_1, u_n_m_2, u_n_m_3, K, rhs);
        end
            
        %( M_div_h2 + K )\(2 * M_div_h2 * u_n - M_div_h2 * u_n_m_1 + rhs)-u_n_p_1
        err(iter) = norm(r);
        if (iter == 20)
            d = d + 1
        end
        if (err(iter) < tol || iter == 10 || norm(J\r) > norm(Delta))
            err_arr(i) = norm(r);
            break;
        end
    end
    %write to output csv
    iterationsRES(i,1) = iter;
    iterationsRES(i,2) = err(iter);
    
    if (iter >3)
        %order(i) = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))));
    end
    %figure;
    %loglog(err);
    %figure;
    %plot(delta);
    order(i) = iter;
    condition = [condition, cond(J)];
    if bdf2
        u_n_m_3 = u_n_m_2;
        u_n_m_2 = u_n_m_1;
        if i==2
            euler=false;
        end
    end
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
    
    y(i,:) = u_n_p_1;
%    if i>=3900
%        disp(u_n_p_1)
%    end
end
toc
d
%order
figure;
plot(t,y(:,1));
title('car position');
%plot(t,y(:,4:11))
figure;
plot(order);
title('newton iterations');
legend()
figure;
plot(err_arr);
title('error');
legend();
%dlmwrite('lookUp11Dof.txt',[y,err_arr],'delimiter',',','precision',9);
dlmwrite('Newton11DOF.txt',iterationsRES,'delimiter',',','precision',12);

function K = get_K()
    global l_long_fl;
    global l_long_fr;
    global l_long_rl;
    global l_long_rr;
    global l_lat_fl;
    global l_lat_fr;
    global l_lat_rl;
    global l_lat_rr;
    global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
    global l1 l2 l3 l4 l5 l6 l7 l8
    
    k1=fnval(k_spline1,l1);
    k2=fnval(k_spline2,l2);
    k3=fnval(k_spline3,l3);
    k4=fnval(k_spline4,l4);
    k5=fnval(k_spline5,l5);
    k6=fnval(k_spline6,l6);
    k7=fnval(k_spline7,l7);
    k8=fnval(k_spline8,l8);
    
    K = [k1+k3+k5+k7, k1*l_lat_fl-k3*l_lat_fr+k5*l_lat_rl-k7*l_lat_rr, -k1*l_long_fl-k3*l_long_fr+k5*l_long_rl+k7*l_long_rr,  -k1, 0, -k3, 0, -k5, 0, -k7, 0;
   0, l_lat_fl*l_lat_fl*k1+l_lat_fr*l_lat_fr*k3+l_lat_rl*l_lat_rl*k5+l_lat_rr*l_lat_rr*k7, -l_long_fl*l_lat_fl*k1+l_lat_fr*l_long_fr*k3+l_long_rl*l_lat_rl*k5-l_long_rr*l_lat_rr*k7, -l_lat_fl*k1, 0, l_lat_fr*k3, 0, -l_lat_rl*k5, 0, l_lat_rr*k7, 0;
   0, 0, l_long_fl*l_long_fl*k1+l_long_fr*l_long_fr*k3+l_long_rl*l_long_rl*k5+l_long_rr*l_long_rr*k7, l_long_fl*k1, 0, l_long_fr*k3, 0, -l_long_rl*k5, 0, -l_long_rr*k7, 0; 
   0, 0, 0, k1+k2, -k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k3+k4, -k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k5+k6, -k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k7+k8, -k8;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k8];
    K=K+K.'-diag(diag(K));
end

function D = get_D()
    global l_long_fl;
    global l_long_fr;
    global l_long_rl;
    global l_long_rr;
    global l_lat_fl;
    global l_lat_fr;
    global l_lat_rl;
    global l_lat_rr;
    global d_spline1 d_spline2 d_spline3 d_spline4 d_spline5 d_spline6 d_spline7 d_spline8;
    global l1 l2 l3 l4 l5 l6 l7 l8
    
    k1=fnval(d_spline1,l1);
    k2=fnval(d_spline2,l2);
    k3=fnval(d_spline3,l3);
    k4=fnval(d_spline4,l4);
    k5=fnval(d_spline5,l5);
    k6=fnval(d_spline6,l6);
    k7=fnval(d_spline7,l7);
    k8=fnval(d_spline8,l8);
    
    K = [k1+k3+k5+k7, k1*l_lat_fl-k3*l_lat_fr+k5*l_lat_rl-k7*l_lat_rr, -k1*l_long_fl-k3*l_long_fr+k5*l_long_rl+k7*l_long_rr,  -k1, 0, -k3, 0, -k5, 0, -k7, 0;
   0, l_lat_fl*l_lat_fl*k1+l_lat_fr*l_lat_fr*k3+l_lat_rl*l_lat_rl*k5+l_lat_rr*l_lat_rr*k7, -l_long_fl*l_lat_fl*k1+l_lat_fr*l_long_fr*k3+l_long_rl*l_lat_rl*k5-l_long_rr*l_lat_rr*k7, -l_lat_fl*k1, 0, l_lat_fr*k3, 0, -l_lat_rl*k5, 0, l_lat_rr*k7, 0;
   0, 0, l_long_fl*l_long_fl*k1+l_long_fr*l_long_fr*k3+l_long_rl*l_long_rl*k5+l_long_rr*l_long_rr*k7, l_long_fl*k1, 0, l_long_fr*k3, 0, -l_long_rl*k5, 0, -l_long_rr*k7, 0; 
   0, 0, 0, k1+k2, -k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k3+k4, -k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k5+k6, -k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k7+k8, -k8;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k8];
    D=K+K.'-diag(diag(K));
end

function dk = getdk_dx()
    dk = 0;
    global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8;
    global d1 d2 d3 d4 d5 d6 d7 d8
    global l1 l2 l3 l4 l5 l6 l7 l8
    d1 = fnval(k_der1, l1);
    d2 = fnval(k_der2, l2);
    d3 = fnval(k_der3, l3);
    d4 = fnval(k_der4, l4);
    d5 = fnval(k_der5, l5);
    d6 = fnval(k_der6, l6);
    d7 = fnval(k_der7, l7);
    d8 = fnval(k_der8, l8);
end

function dd = getdd_dx()
    dk = 0;
    global d_der1 d_der2 d_der3 d_der4 d_der5 d_der6 d_der7 d_der8;
    global d1 d2 d3 d4 d5 d6 d7 d8
    global l1 l2 l3 l4 l5 l6 l7 l8
    d1 = fnval(d_der1, l1);
    d2 = fnval(d_der2, l2);
    d3 = fnval(d_der3, l3);
    d4 = fnval(d_der4, l4);
    d5 = fnval(d_der5, l5);
    d6 = fnval(d_der6, l6);
    d7 = fnval(d_der7, l7);
    d8 = fnval(d_der8, l8);
end