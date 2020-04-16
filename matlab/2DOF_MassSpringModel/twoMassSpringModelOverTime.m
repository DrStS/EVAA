clear; clc; close all;

%%
a = 1;
b = 2;
c = 3;
f1 = 1;
f2 = 0;
rhs = [f1;f2];
L1 = 1;
L2 = 1;

% parameters
tol = 1e-12;
x_init = [0; 0];
err = [];

num_iter = 1000;
delta_t = 1e-1; 
t = 0:delta_t:num_iter*delta_t;

m_1 = 1;
m_2 = 1;
M = diag([m_1,m_2]);
M_div_h2 = M/(delta_t^2);

%%
% evaluation density
% evaluation density
dl = 0.01;

% calc min and max length to evaluate (just for now)
l_min = -1;
l_max = 3;
% grid size
k_d = round((l_max - l_min)/dl);
% grid value allocation
k_grid = zeros(k_d+1,1);
% k response function
k = @(l)(c*l*l + b*l + a);
% store x values
x = zeros(k_d+1,1);
% fill in grid values with chebyshev
% for i = 0:k_d-1
%     x(i+1) = (1+cos((2*i+1)/(2*k_d)*pi))/2*(l_max-l_min)+l_min;
%     k_grid(i+1)= k(x(i+1));
% end
for i = 0:k_d
    x(i+1) = i*dl+l_min;
    k_grid(i+1)= k(x(i+1));
end
% grid spline
k_spline = spline(x,k_grid);
% grid derivative spline
k_der = fnder(k_spline,1);

%% init declare variables

syms l1 u1 u2 l2 ;
l1 = u1 - u2 + L1;
l2 = u2 + L2;

%% Internal forces
syms p1 p2 k1 k2;
p1 = k1*u1 - k1*u2;
p2 = -k1*u1+(k1+k2)*u2;

%% Jacobian
syms K k11 k12 k21 k22 K_symb d1 d2
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

K_symb = [k11, k12; 
    k21, k22];

dKdx_x_symb = [d1*(u1 - u2)	d1*(-u1 + u2);...
    d1*(-u1 + u2)	d1*(u1 - u2) + d2*u2];

f_newton = @(y_curr,y1,y2,K)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);

 
%% Newton iteration
%init variables for symb calc
u1 = 0;
u2 = 0;
u = [u1; u2];
Delta = [u1; u2];

y = zeros(2, length(t));
err_arr = zeros(1,length(t));

%init sol vecotr
y(:,1) = u - 20 * delta_t;
y(:,2) = u;
iterations = [];

% inital values
l_1 = eval(l1);
l_2 = eval(l2);
k1 = fnval(k_spline,l_1);
k2 = fnval(k_spline,l_2);

for i = 2: length(t)-1
    i
    y(:,i+1) = y(:,i);
    K = eval(K_symb);
    r = f_newton(y(:,i+1), y(:,i), y(:,i-1),K);
    
    iter = 0;
    err = [];
    while 1
       iter = iter + 1;
       d1 = fnval(k_der, l_1);
       d2 = fnval(k_der, l_2);
       J = eval(dKdx_x_symb);
       J = M_div_h2 + K + J;
       Delta = -J\r;
       y(:,i+1) = Delta + y(:,i+1);
       % update values
       u1 = y(1,i+1);
       u2 = y(2,i+1);
       l_1 = eval(l1);
       l_2 = eval(l2);
       k1 = fnval(k_spline,l_1);
       k2 = fnval(k_spline,l_2);
       K = eval(K_symb);
       r = f_newton(y(:,i+1), y(:,i), y(:,i-1),K);
       err = [err, norm(r)];
       if (err(iter) < tol)
           err_arr(i) = norm(r);
            break;
       end
    end
    figure;
    loglog(err);
    if iter > 3
        order = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))))
    end
    iterations = [iterations, iter];
end
figure;
title('newton iterations');
plot(iterations);
%plot(t,y);
figure;
title('error');
semilogy(err_arr);