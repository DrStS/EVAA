clear; clc; close all;
%% make lookup
a = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3] ;
b = 1.59;
c = 1.2;

%%
global k_grid size_grid l_min l_max dl;
% calc min and max length to evaluate (just for now)
l_min = 0.05;
l_max = 0.8;
% grid size
size_grid = 1000;
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
%% symbolic derivative
syms l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr
k = sym('k',[8 1]); 
syms k_body_fl k_tyre_fl k_body_fr k_tyre_fr k_body_rl k_tyre_rl k_body_rr k_tyre_rr;
syms K M_div_h2

K = [k(1)+k(3)+k(5)+k(7), k(1)*l_lat_fl-k(3)*l_lat_fr+k(5)*l_lat_rl-k(7)*l_lat_rr, -k(1)*l_long_fl-k(3)*l_long_fr+k(5)*l_long_rl+k(7)*l_long_rr,  -k(1), 0, -k(3), 0, -k(5), 0, -k(7), 0;
   0, l_lat_fl*l_lat_fl*k(1)+l_lat_fr*l_lat_fr*k(3)+l_lat_rl*l_lat_rl*k(5)+l_lat_rr*l_lat_rr*k(7), -l_long_fl*l_lat_fl*k(1)+l_lat_fr*l_long_fr*k(3)+l_long_rl*l_lat_rl*k(5)-l_long_rr*l_lat_rr*k(7), -l_lat_fl*k(1), 0, l_lat_fr*k(3), 0, -l_lat_rl*k(5), 0, l_lat_rr*k(7), 0;
   0, 0, l_long_fl*l_long_fl*k(1)+l_long_fr*l_long_fr*k(3)+l_long_rl*l_long_rl*k(5)+l_long_rr*l_long_rr*k(7), l_long_fl*k(1), 0, l_long_fr*k(3), 0, -l_long_rl*k(5), 0, -l_long_rr*k(7), 0; 
   0, 0, 0, k(1)+k(2), -k(2), 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k(2), 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k(3)+k(4), -k(4), 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k(4), 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k(5)+k(6), -k(6), 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k(6), 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k(7)+k(8), -k(8);
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k(8)];

K=simplify(K+K.'-diag(diag(K)));

dK1_dk = [diff(K(:,1), k(1)), diff(K(:,1), k(2)), diff(K(:,1), k(3)), diff(K(:,1), k(4)), diff(K(:,1), k(5)), diff(K(:,1), k(6)), diff(K(:,1), k(7)), diff(K(:,1), k(8))];
dK2_dk = [diff(K(:,2), k(1)), diff(K(:,2), k(2)), diff(K(:,2), k(3)), diff(K(:,2), k(4)), diff(K(:,2), k(5)), diff(K(:,2), k(6)), diff(K(:,2), k(7)), diff(K(:,2), k(8))];
dK3_dk = [diff(K(:,3), k(1)), diff(K(:,3), k(2)), diff(K(:,3), k(3)), diff(K(:,3), k(4)), diff(K(:,3), k(5)), diff(K(:,3), k(6)), diff(K(:,3), k(7)), diff(K(:,3), k(8))];
dK4_dk = [diff(K(:,4), k(1)), diff(K(:,4), k(2)), diff(K(:,4), k(3)), diff(K(:,4), k(4)), diff(K(:,4), k(5)), diff(K(:,4), k(6)), diff(K(:,4), k(7)), diff(K(:,4), k(8))];
dK5_dk = [diff(K(:,5), k(1)), diff(K(:,5), k(2)), diff(K(:,5), k(3)), diff(K(:,5), k(4)), diff(K(:,5), k(5)), diff(K(:,5), k(6)), diff(K(:,5), k(7)), diff(K(:,5), k(8))];
dK6_dk = [diff(K(:,6), k(1)), diff(K(:,6), k(2)), diff(K(:,6), k(3)), diff(K(:,6), k(4)), diff(K(:,6), k(5)), diff(K(:,6), k(6)), diff(K(:,6), k(7)), diff(K(:,6), k(8))];
dK7_dk = [diff(K(:,7), k(1)), diff(K(:,7), k(2)), diff(K(:,7), k(3)), diff(K(:,7), k(4)), diff(K(:,7), k(5)), diff(K(:,7), k(6)), diff(K(:,7), k(7)), diff(K(:,7), k(8))];
dK8_dk = [diff(K(:,8), k(1)), diff(K(:,8), k(2)), diff(K(:,8), k(3)), diff(K(:,8), k(4)), diff(K(:,8), k(5)), diff(K(:,8), k(6)), diff(K(:,8), k(7)), diff(K(:,8), k(8))];
dK9_dk = [diff(K(:,9), k(1)), diff(K(:,9), k(2)), diff(K(:,9), k(3)), diff(K(:,9), k(4)), diff(K(:,9), k(5)), diff(K(:,9), k(6)), diff(K(:,9), k(7)), diff(K(:,9), k(8))];
dK10_dk = [diff(K(:,10), k(1)), diff(K(:,10), k(2)), diff(K(:,10), k(3)), diff(K(:,10), k(4)), diff(K(:,10), k(5)), diff(K(:,10), k(6)), diff(K(:,10), k(7)), diff(K(:,10), k(8))];
dK11_dk = [diff(K(:,11), k(1)), diff(K(:,11), k(2)), diff(K(:,11), k(3)), diff(K(:,11), k(4)), diff(K(:,11), k(5)), diff(K(:,11), k(6)), diff(K(:,11), k(7)), diff(K(:,11), k(8))];
dKcols_dk = [dK1_dk, dK2_dk, dK3_dk, dK4_dk, dK5_dk, dK6_dk, dK7_dk, dK8_dk, dK9_dk, dK10_dk, dK11_dk];



%% parameters
% time
num_iter = 1000;
delta_t = 1e-3; 
t = 0:delta_t:(num_iter-1)*delta_t;


tol = 1e-12;
y = zeros(11, length(t));
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
g = 0;
u_n = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

y(:,1) = u_n;
u_n_m_1=u_n;
u_n_p_1=u_n;

rhs =[1.1e3; zeros(10,1)];
M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
M_div_h2 = M / (delta_t * delta_t);

%% time steps

f_newton = @(y_curr,y1,y2,K)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);

dKcols_dk = eval(dKcols_dk); % evaluate it numerically

d = 0;
for i = 1: length(t)-1
    K = get_K(u_n);
    r = f_newton(u_n_p_1, u_n, u_n_m_1, K);
    
    iter = 0;
    err = [];
    while 1
        iter = iter + 1;
        temp = [];
        dk_dx = getdk_dx(u_n_p_1);
        for ii = 1 : 11 % 11 columns in K
            temp = [temp; dk_dx * u_n_p_1(ii)];
        end
        J = M_div_h2 + K + dKcols_dk * temp;
        
        Delta = -J\r;
%         norm(J)
        u_n_p_1 = Delta + u_n_p_1; 
        
        % update values
        K = get_K(u_n_p_1);
        r = f_newton(u_n_p_1, u_n, u_n_m_1, K);
        err(iter) = norm(r);
        if (iter == 10)
            d = d + 1;
        end
        if (err(iter) < tol || iter == 10)
            err_arr(i) = norm(r);
            break;
        end
    end
    
    if (iter >3)
        order(i) = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))));
    end
    
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
    y(:,i+1) = u_n_p_1;
end
d
order
%plot(t,y(1,:));
plot(t,y(4:11,:))
%plot(err_arr);
legend();

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
    global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
    
    %global X k_grid size_grid;
    %k_body_fl=interp1(X,k_grid(1:size_grid),0.5 - upper_fl + car);
    %k_tyre_fl=interp1(X,k_grid(size_grid+1:2*size_grid),0.5 + upper_fl - lower_fl);
    %k_body_fr=interp1(X,k_grid(2*size_grid+1:3*size_grid),0.5 - upper_fr + car);
    %k_tyre_fr=interp1(X,k_grid(3*size_grid+1:4*size_grid),0.5 + upper_fr - lower_fr);
    %k_body_rl=interp1(X,k_grid(4*size_grid+1:5*size_grid),0.5 - upper_rl + car);
    %k_tyre_rl=interp1(X,k_grid(5*size_grid+1:6*size_grid),0.5 + upper_rl -lower_rl);
    %k_body_rr=interp1(X,k_grid(6*size_grid+1:7*size_grid),0.5 - upper_rr + car);
    %k_tyre_rr=interp1(X,k_grid(7*size_grid+1:8*size_grid),0.5 + upper_rr - lower_rr);
    
    k_body_fl=fnval(k_spline1,0.5 - upper_fl + car);
    k_tyre_fl=fnval(k_spline2,0.5 + upper_fl - lower_fl);
    k_body_fr=fnval(k_spline3,0.5 - upper_fr + car);
    k_tyre_fr=fnval(k_spline4,0.5 + upper_fr - lower_fr);
    k_body_rl=fnval(k_spline5,0.5 - upper_rl + car);
    k_tyre_rl=fnval(k_spline6,0.5 + upper_rl -lower_rl);
    k_body_rr=fnval(k_spline7,0.5 - upper_rr + car);
    k_tyre_rr=fnval(k_spline8,0.5 + upper_rr - lower_rr);
    
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
    K=K+K.'-diag(diag(K));
end

function dk = getdk_dx(x) % 8x11
    dk = zeros(8,11);
    
    car = x(1);
    upper_fl = x(4);
    lower_fl = x(5);
    upper_fr = x(6);
    lower_fr = x(7);
    upper_rl = x(8);
    lower_rl = x(9);
    upper_rr = x(10);
    lower_rr = x(11);
    
    %dk1_dl = deriv_dk(x(4)-x(5), 1);
    %dk2_dl = deriv_dk(x(5), 2);
    %dk3_dl = deriv_dk(x(6)-x(7), 3);
    %dk4_dl = deriv_dk(x(7), 4);
    %dk5_dl = deriv_dk(x(8)-x(9), 5);
    %dk6_dl = deriv_dk(x(9), 6);
    %dk7_dl = deriv_dk(x(10)-x(11), 7);
    %dk8_dl = deriv_dk(x(11), 8);
    
    global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 
    
    dk1_dl = fnval(k_der1, 0.5 - upper_fl + car);
    dk2_dl = fnval(k_der2, 0.5 + upper_fl - lower_fl);
    dk3_dl = fnval(k_der3, 0.5 - upper_fr + car);
    dk4_dl = fnval(k_der4, 0.5 + upper_fr - lower_fr);
    dk5_dl = fnval(k_der5, 0.5 - upper_rl + car);
    dk6_dl = fnval(k_der6, 0.5 + upper_rl -lower_rl);
    dk7_dl = fnval(k_der7, 0.5 - upper_rr + car);
    dk8_dl = fnval(k_der8, 0.5 + upper_rr - lower_rr);
    
    dk(1,1) = dk1_dl;
    dk(1,4) = -dk1_dl;
    dk(2,4) = dk2_dl;
    dk(2,5) = -dk2_dl;
    dk(3,1) = dk3_dl;
    dk(3,6) = -dk3_dl;
    dk(4,6) = dk4_dl;
    dk(4,7) = -dk4_dl;
    dk(5,1) = dk5_dl;
    dk(5,8) = -dk5_dl;
    dk(6,8) = dk6_dl;
    dk(6,9) = -dk6_dl;
    dk(7,1) = dk7_dl;
    dk(7,10) = -dk7_dl;
    dk(8,10) = dk8_dl;
    dk(8,11) = -dk8_dl;
end
