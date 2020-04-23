clear; clc; close all;
%% make lookup
a = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3] ;
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
global k_grid size_grid l_min l_max dl;
% calc min and max length to evaluate (just for now)
l_min = 0.02;
l_max = 1;
% grid size
size_grid = 1001;
% grid value allocation
global X;
X = zeros(size_grid, 1);
k_grid = zeros(8*size_grid,1);
dl = (l_max-l_min)/(size_grid-1);
% k response function
k = @(l,a)(c*l*l + b*l + a);
% fill in grid values
for i = 0:size_grid-1
    X(i+1) = l_min+i*dl;
end
for i = 1:8
    for j = 1:size_grid
        k_grid((i-1)*size_grid + j)= k(X(j), a(i));
    end
end
global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 
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
num_iter = 4;
delta_t = 1e-1; 
t = 0:delta_t:(num_iter-1)*delta_t;

tol = 1e-10;
y = zeros(length(t),7);
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
I_body_xx=640;
I_body_yy=4800;
mass_wheel_fl=67.5;
mass_tyre_fl=30;
mass_wheel_fr=67.5;
mass_tyre_fr=30;
mass_wheel_rl=67.5;
mass_tyre_rl=30;
mass_wheel_rr=67.5;
mass_tyre_rr=30;
u_n = [0; 0; 0; 0; 0; 0; 0];
%% init values
x1 = u_n(1);
x2 = u_n(2);
x3 = u_n(3);
x4 = u_n(4);
x5 = 0;
x6 = u_n(5);
x7 = 0;
x8 = u_n(6);
x9 = 0;
x10 = u_n(7);
x11 = 0;
l1 = eval(l1_sym);
l2 = eval(l2_sym);
l3 = eval(l3_sym);
l4 = eval(l4_sym);
l5 = eval(l5_sym);
l6 = eval(l6_sym);
l7 = eval(l7_sym);
l8 = eval(l8_sym);

u_n_m_1 = u_n;
u_n_m_1(1)=u_n(1) - 0 * delta_t;
u_n_p_1=u_n;

idx = [5 7 9 11];
rhs =-[mass_Body; 0; 0; mass_wheel_fl; mass_wheel_fr; mass_wheel_rl; mass_wheel_rr];
%rhs =-[mass_Body; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
%rhs = 0;
M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_wheel_fr, mass_wheel_rl, mass_wheel_rr]);
M_div_h2 = M / (delta_t * delta_t);

%% time steps

f_newton = @(y_curr,y1,y2,K)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);
%dKcols_dk = eval(dKcols_dk); % evaluate it numerically
condition = [];
d = 0; grad_cond = 0;
for i = 1: length(t)
    K = get_K();
    K(:,idx) = [];
    K(idx,:) = [];
    new_r = f_newton(u_n_p_1, u_n, u_n_m_1, K);
    i
    iter = 0;
    err = [];
    delta = [];
    init_err = new_r;
    while 1
        r = new_r;
        iter = iter + 1;
        temp = [];
        getdk_dx();
        dKdx_x = eval(dKdx_x_symb);
        dKdx_x(:,idx) = [];
        dKdx_x(idx,:) = [];
        J = M_div_h2 + K + dKdx_x;
        Delta = -J\r;
        %delta(iter) = norm(Delta);
        u_n_p_1 = Delta + u_n_p_1; 
        %u_n_p_1(idx) = 0;
        
        % update values
        x1 = u_n_p_1(1);
        x2 = u_n_p_1(2);
        x3 = u_n_p_1(3);
        x4 = u_n_p_1(4);
        x6 = u_n_p_1(5);
        x8 = u_n_p_1(6);
        x10 = u_n_p_1(7);
        l1 = eval(l1_sym);
        l2 = eval(l2_sym);
        l3 = eval(l3_sym);
        l4 = eval(l4_sym);
        l5 = eval(l5_sym);
        l6 = eval(l6_sym);
        l7 = eval(l7_sym);
        l8 = eval(l8_sym);
        K = get_K();
        K(:,idx) = [];
        K(idx,:) = [];
        new_r = f_newton(u_n_p_1, u_n, u_n_m_1, K);
        %( M_div_h2 + K )\(2 * M_div_h2 * u_n - M_div_h2 * u_n_m_1 + rhs)-u_n_p_1
        err(iter) = norm(new_r);
        if (iter == 20)
            d = d + 1;
        end
        if (err(iter) < tol || (norm(J\new_r) > norm(Delta)) || iter == 20)
            if((norm(J\new_r) > norm(Delta)))
              grad_cond = grad_cond + 1;
            end
            err_arr(i) = norm(new_r);
            break;
        end
    end
    
    if (iter >3)
        %order(i) = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))));
    end
    %figure;
    %loglog(err);
    %figure;
    %plot(delta);
    order(i) = iter;
    condition = [condition, cond(J)];
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
    y(i,:) = u_n_p_1;
end
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
%csvwrite('lookUp11DofReduced.txt',[y,err_arr]);
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
