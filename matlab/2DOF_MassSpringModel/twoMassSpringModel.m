clear; clc; close all;

%% parameters
m1 = 1; m2 = 2;
k1_init = 350; k2_init = 100;
g = 9.81;
l1 = 1;
l2 = 1;

tol = 1e-8;


%% generate model (k1,k2)=f(x1,x2)
interpolate_k = @(x)([2, 1; 1, 5]*[(l1-x(1))^2; (l2-x(2)+x(1))^2]+[k1_init; k2_init]);
                  
%% declare aux_vals
aux_vals = struct('m1', m1, ...
                  'm2', m2, ...
                  'l1', l1, ...
                  'l2', l2, ...
                  'g',  g, ...
                  'interpolate_k', interpolate_k, ...
                  'tol', tol);

x_prev = [l1; l1+l2]; 
x_curr = x_prev;
              
num_iter = 1000;
delta_t = 1e-3; 
t = 0:delta_t:num_iter*delta_t;


%% equilibrium

[y2,k2] = newton_equi(x_curr, aux_vals);
d1 = (l1-y2(1));
d2 = (l2-y2(2)+y2(1));
abs(d1*k2(1) - (m1+m2)*g) < tol
abs(d2*k2(2) - m2*g) < tol

%% simulation

[t,y] = BW_2DOF_scheme(t, x_prev, x_curr, aux_vals);

%% plots
figure;
plot(t,y);
hold on;
y_equi = y2 * ones(1, length(y(1,:)));
plot(t,y_equi);
grid on;
legend('x1','x2','x1 equi', 'x2 equi');


