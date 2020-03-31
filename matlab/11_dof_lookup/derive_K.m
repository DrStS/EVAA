clear; clc; close all;
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

lookup_struct = load('lookup_table.mat');
lookup_table = lookup_struct.big_lookup_table;
global k_grid size_grid l_min l_max dl;
size_grid = length(lookup_table);
l_min = lookup_table(1,1);
l_max = lookup_table(size_grid,1);
dl = (l_max - l_min)/(size_grid-1);
k_grid = lookup_table(:,2:9);

%% parameters
% time
num_iter = 1000;
delta_t = 1e-5; 
t = 0:delta_t:num_iter*delta_t;


tol = 1e-12;
y = zeros(11, length(t));
order = zeros(length(t),1);
err_arr = zeros(length(t),1);

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
g = 1;
u_n = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
u_n = rand(11,1);

y(:,1) = u_n;
u_n_m_1=u_n;
u_n_p_1=u_n;

rhs = [mass_Body; 0; 0; mass_wheel_fl; mass_tyre_fl; mass_wheel_fr; mass_tyre_fr; mass_wheel_rl; mass_tyre_rl; mass_wheel_rr; mass_tyre_rr]  * g;
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
    upper_fl = x(4);
    lower_fl = x(5);
    upper_fr = x(6);
    lower_fr = x(7);
    upper_rl = x(8);
    lower_rl = x(9);
    upper_rr = x(10);
    lower_rr = x(11);
    
    k_body_fl=interpolate(upper_fl - lower_fl, 1);
    k_tyre_fl=interpolate(lower_fl, 2);
    k_body_fr=interpolate(upper_fr-lower_fr, 3);
    k_tyre_fr=interpolate(lower_fr, 4);
    k_body_rl=interpolate(upper_rl-lower_rl, 5);
    k_tyre_rl=interpolate(lower_rl, 6);
    k_body_rr=interpolate(upper_rr-lower_rr, 7);
    k_tyre_rr=interpolate(lower_rr, 8);
    
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
    
    dk1_dl = deriv_dk(x(4)-x(5), 1);
    dk2_dl = deriv_dk(x(5), 2);
    dk3_dl = deriv_dk(x(6)-x(7), 3);
    dk4_dl = deriv_dk(x(7), 4);
    dk5_dl = deriv_dk(x(8)-x(9), 5);
    dk6_dl = deriv_dk(x(9), 6);
    dk7_dl = deriv_dk(x(10)-x(11), 7);
    dk8_dl = deriv_dk(x(11), 8);
    
    dk(1,4) = dk1_dl;
    dk(1,5) = -dk1_dl;
    dk(2,5) = dk2_dl;
    dk(3,6) = dk3_dl;
    dk(3,7) = -dk3_dl;
    dk(4,7) = dk4_dl;
    dk(5,8) = dk5_dl;
    dk(5,9) = -dk5_dl;
    dk(6,9) = dk6_dl;
    dk(7,10) = dk7_dl;
    dk(7,11) = -dk7_dl;
    dk(8,11) = dk8_dl;
end

% interpolates k_grid to a certain value l
function k = interpolate(l, s)
    global k_grid l_min l_max dl;
    if (l < l_min)
        l = l_min;
    end
    if (l > l_max) 
       l = l_max; 
    end
    % matlab need a +1 
    % we dont need it in c++ tho
    idx = (l-l_min)/dl + 1;
    if (floor(idx) == idx)
        k = k_grid(idx,s);
    else
        i_low = floor(idx);
        low = k_grid(i_low, s);
        top = k_grid(ceil(idx), s);
        pr = (l-(l_min+dl*(i_low-1)))/dl; % percentage from left to right
        k = (1-pr)*low+pr*top;
    end
end

% calculates the derivative at a certain position of the grid
function dk = deriv_dk(l, s)
    global k_grid l_min l_max dl;
    if l <= l_min
       l = l_min + 1e-12; 
    end
    if l >= l_max
       l = l_max - 1e-12; 
    end
    idx = (l-l_min)/dl + 1;
    if floor(idx) == idx
        low = k_grid(idx-1, s);
        top = k_grid(idx+1, s);
        dk = (top-low)/(2*dl);
    else
        low = k_grid(floor(idx), s);
        top = k_grid(ceil(idx), s);
        dk = (top-low)/dl;
    end
end