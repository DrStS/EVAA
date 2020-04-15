clear; clc; close all;

%%
a = 1;
b = 2;
c = 3;
f1 = 0;
f2 = 1;
L1 = 1;
L2 = 1;

% parameters
tol = 1e-12;
x_init = [0; 0];
iter = 0;
err = [];

%% interpolation
% evaluation density
dl = 0.01;

% calc min and max length to evaluate (just for now)
l_min = 0.05;
l_max = 1.7;
% grid size
k_d = round((l_max - l_min)/dl);
% grid value allocation
k_grid = zeros(k_d+1,1);
% k response function
k = @(l)(c*l*l + b*l + a);
% store x values
x = zeros(k_d+1,1);
% fill in grid values with chebyshev
%for i = 0:k_d-1
%    x(i+1) = (1+cos((2*i+1)/(2*k_d)*pi))/2*(l_max-l_min)+l_min;
%    k_grid(i+1)= k(x(i+1));
%end
for i = 0:k_d
    x(i+1) = i*dl+l_min;
    k_grid(i+1)= k(x(i+1));
end
%plot(x,k_grid);
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
p1 = k1*u1 - k2*u2;
p2 = -k1*u1+(k1+k2)*u2;

%% Jacobian
syms K k11 k12 k21 k22 K_symb d1 d2
% actuall derivative but sym cant be used as we get part of derivative out
% of the spline
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

r_symb = [p1 - f1; 
     p2 - f2];

dKdx_x_symb = [d1*(u1 - u2)	d1*(-u1 + u2);...
    d1*(-u1 + u2)	d1*(u1 - u2) + d2*u2];

 
%% Newton iteration
u1 = 0;
u2 = 0;
u = [u1; u2];
Delta = [u1; u2];

% inital values
l_1 = eval(l1);
l_2 = eval(l2);
k1 = fnval(k_spline,l_1);
k2 = fnval(k_spline,l_2);
K = eval(K_symb);
r = eval(r_symb);
while 1
   iter = iter + 1;
   d1 = fnval(k_der, l_1);
   d2 = fnval(k_der, l_2);
   J = eval(dKdx_x_symb);
   J = K + J;
   Delta = -J\r;
   u = Delta + u;
   % update values
   u1 = u(1);
   u2 = u(2);
   l_1 = eval(l1);
   l_2 = eval(l2);
   k1 = fnval(k_spline,l_1);
   k2 = fnval(k_spline,l_2);
   K = eval(K_symb);
   r = eval(r_symb);
   err = [err, norm(r)];
   if (err(length(err)) < tol)
        break;
   end
end
iter
order = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))))
disp(u)
disp(err(length(err)))
ref = [ 0.2788; 0.1394]
figure;
semilogy(err)
title('error');