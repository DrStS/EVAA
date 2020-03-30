clear; clc; close all;

%%
a = 1;
b = 2;
c = 3;
f1 = -1;
f2 = -1;
rhs = [f1;f2];
L1 = 1;
L2 = 1;

% parameters
tol = 1e-12;
x_init = [0; 0];
err = [];

num_iter = 1000;
delta_t = 5e-4; 
t = 0:delta_t:num_iter*delta_t;

m_1 = 1;
m_2 = 1;
M = diag([m_1,m_2]);
M_div_h2 = M/delta_t^2;

%%
% evaluation density
dl = 0.1;
% calc min and max length to evaluate (just for now)
l_min = min(0.3 * L1,0.3 * L2);
l_max = max(1.7 * L1,1.7 * L2);
% grid size
k_d = round((l_max - l_min)/dl);
% grid value allocation
k_grid = zeros(k_d+1,1);
% k response function
k = @(l)(c*l*l + b*l + a);
% fill in grid values
for i = 0:k_d
        k_grid(i+1)= k(l_min+i*dl);
end

%% init declare variables

syms l1 u1 u2 l2 ;
l1 = u1 - u2 + L1;
l2 = u2 + L2;

%% Internal forces
syms p1 p2 k1 k2;
p1 = k1*u1 - k2*u2;
p2 = -k1*u1+(k1+k2)*u2;

%% Jacobian
syms K k11 k12 k21 k22 K_symb
% actuall derivative but sym cant be used as we get part of derivative out
% my function
%k11 = diff(p1, u1)+u1*deriv_dk(l1);
%k12 = diff(p1, u2)-u1*deriv_dk(l1)-u2*deriv_dk(l2);
%k21 = diff(p2, u1)+(u2-u1)*deriv_dk(l1);
%k22 = diff(p2, u2)-(u2-u1)*deriv_dk(l1)+u2*deriv_dk(l2);
k11 = diff(p1, u1);
k12 = diff(p1, u2);
k21 = diff(p2, u1);
k22 = diff(p2, u2);

K_symb_1 = [k11, k12; 
    k21, k22];

f_newton = @(y_curr,y1,y2,K)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);

 
%% Newton iteration
u1 = 0;
u2 = 0;
u = [u1; u2];
Delta = [u1; u2];

y = zeros(2, length(t));
err = zeros(1,length(t));
y(:,1) = u;
y(:,2) = u;

% inital values
l_1 = eval(l1);
l_2 = eval(l2);
k1 = interpolate(l_1,l_min,dl,k_grid);
k2 = interpolate(l_2,l_min,dl,k_grid);

for i = 2: length(t)-1
    y(:,i+1) = y(:,i);
    K = eval(K_symb_1)+[y(1,i+1),-y(1,i+1);y(2,i+1)-y(1,i+1),y(1,i+1)-y(2,i+1)]*deriv_dk(l_1,l_min,dl,k_grid)+ [0 -y(2,i+1); 0 y(2,i+1)]*deriv_dk(l_2,l_min,dl,k_grid);
    r = f_newton(y(:,i+1), y(:,i), y(:,i-1),K);
    
    iter = 0;
    while 1
        iter = iter + 1;
        J = M_div_h2 + K;
        Delta = -J\r;
        y(:,i+1) = Delta + y(:,i+1);
        % update values
        l_1 = eval(l1);
        l_2 = eval(l2);
        k1 = interpolate(l_1,l_min,dl,k_grid);
        k2 = interpolate(l_2,l_min,dl,k_grid);
        K = eval(K_symb_1)+[y(1,i+1),-y(1,i+1);y(2,i+1)-y(1,i+1),y(1,i+1)-y(2,i+1)]*deriv_dk(l_1,l_min,dl,k_grid)+ [0 -y(2,i+1); 0 y(2,i+1)]*deriv_dk(l_2,l_min,dl,k_grid);
        r = f_newton(y(:,i+1), y(:,i), y(:,i-1),K);
        err(i) = norm(r);
        if (err(i) < tol || iter == 10)
            break;
        end
    end
end

plot(t,y);
plot(err);

% interpolates k_grid to a certain value l
function k = interpolate(l,l_min,dl,k_grid)
    i_low = floor((l-l_min)/dl);
    low = k_grid(i_low);
    top = k_grid(ceil((l-l_min)/dl));

    pr = (l-l_min-dl*i_low)/dl; % percentage from left to right
    k = (1-pr)*low+pr*top;
end

% calculates the derivative at a certain position of the grid
function dk = deriv_dk(l,l_min,dl,k_grid)
    low = k_grid(floor((l-l_min)/dl));
    top = k_grid(ceil((l-l_min)/dl));
    dk = (top-low)/dl;
end