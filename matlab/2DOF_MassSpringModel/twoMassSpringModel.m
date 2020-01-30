clear; clc; close all;

%% parameters
m1 = 1; m2 = 1*2;
k1_init = 350; k2_init = 100;
g = 9.81;
l1 = 1;
l2 = 1;

tol = 1e-8;


%% generate model (k1,k2)=f(x1,x2)
quadratic_part = [2, 1; 1, 5];
 lin_part = [-2 -3; -1, -5];
%lin_part = [0, 0; 0, 0];

interpolate_k = @(dx)(quadratic_part*[dx(1)^2; dx(2)^2]+lin_part*[dx(1);dx(2)]+[k1_init; k2_init]);
% interpolate_k = @(dx)([2, 1; 1, 5]*[(l1-x(1))^2; (l2-x(2)+x(1))^2]+[k1_init; k2_init]);
dk_d_dx = @(dx)(2*quadratic_part*diag(dx) + lin_part);

f_rhs = @(k)([-m1*g + k(1)*l1 - k(2)*l2; ...
              -m2*g + k(2)*l2]);
df_rhs_dk = [l1, -l2; 0, l2];  
        
f_K = @(k)([k(1)+k(2), -k(2); 
           -k(2), k(2)]);
dK_dk1 = [1, 0; 0, 0];
dK_dk2 = [1, -1; -1, 1];
       
get_dx = @(x)([l1-x(1); 
               l2+x(1)-x(2)]);
ddx_dx = [-1, 0; ...
           1, -1];           
           
%% declare aux_vals
aux_vals = struct('m1', m1, ...
                  'm2', m2, ...
                  'l1', l1, ...
                  'l2', l2, ...
                  'g',  g, ...
                  'interpolate_k', interpolate_k, ...
                  'f_K', f_K, ...
                  'f_rhs', f_rhs, ...
                  'dK_dk1', dK_dk1, ...
                  'dK_dk2', dK_dk2, ...
                  'dk_d_dx', dk_d_dx, ...
                  'ddx_dx', ddx_dx, ...
                  'get_dx', get_dx, ...
                  'df_rhs_dk', df_rhs_dk, ...
                  'tol', tol);

x_prev = [l1; l1+l2]; 
x_curr = x_prev;
              
num_iter = 1000;
delta_t = 1e-3; 
t = 0:delta_t:num_iter*delta_t;


%% equilibrium

[y2,k2] = Copy_of_newton_equi(x_curr, aux_vals);
d1 = (l1-y2(1));
d2 = (l2-y2(2)+y2(1));
disp('Equilibrium solution');
abs(d1*k2(1) - (m1+m2)*g) < tol
abs(d2*k2(2) - m2*g) < tol

%% simulation

[t,y] = Copy_of_BW_2DOF_scheme(t, x_prev, x_curr, aux_vals);

%% plots
figure;
plot(t,y);
hold on;
y_equi = y2 * ones(1, length(y(1,:)));
plot(t,y_equi);
grid on;
legend('x1','x2','x1 equi', 'x2 equi');


