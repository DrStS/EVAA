clear; clc; close all;

%% parameters
m1 = .1; m2 = .2;
k1_init = 10; k2_init = 5;
g = 9.81;

tol = 1e-10;

x1 = 0;
x2 = 0;
x_init = [0; 0];
              
num_iter = 1000;
delta_t = 5e-4; 
t = 0:delta_t:num_iter*delta_t;

%% jacobi without partial derivatives 
Q11 = 2;
Q12 = 1;
Q21 = 1;
Q22 = 5;

Jacobi_equi = @(x)([ 3*Q11*x(1)^2 - 2*Q21*x(1)*x(2) + Q12*x(2)^2 + k1_init, - Q21*x(1)^2 + 2*Q12*x(1)*x(2) - 3*Q22*x(2)^2 - k2_init; ...
           2*Q21*x(1)*x(2),                 Q21*x(1)^2 + 3*Q22*x(2)^2 + k2_init]);
       
Jacobi_time = @(x)([ m1/delta_t^2 + 3*Q11*x(1)^2 - 2*Q21*x(1)*x(2) + Q12*x(2)^2 + k1_init, - Q21*x(1)^2 + 2*Q12*x(1)*x(2) - 3*Q22*x(2)^2 - k2_init; ...
           2*Q21*x(1)*x(2),                 m2/delta_t^2 + Q21*x(1)^2 + 3*Q22*x(2)^2 + k2_init]);
       
%% generate model (k1,k2)=f(x1,x2)

rhs = [-m1*g; ...
         -m2*g ];
     
quadratic_part = [2, 1; 1, 5];
lin_part = [-2 -3; -1, -5];
lin_part = [0, 0; 0, 0];
%quadratic_part = [0, 0; 0, 0];

f_k = @(x)(quadratic_part*[x(1)^2; x(2)^2]+lin_part*[x(1);x(2)]+[k1_init; k2_init]);
f_dk_dx = @(x)(2*quadratic_part*diag(x) + lin_part);
        
f_K = @(k)([k(1), -k(2); 
           0, k(2)]);
dK_dk1 = [1, 0; 0, 0];
dK_dk2 = [0, -1; 0, 1];

f_tens_vec = @(A1,A2,v)([A1(1,1)*v(1)+A2(1,1)*v(2), A1(1,2)*v(1)+A2(1,2)*v(2);...
    A1(2,1)*v(1)+A2(2,1)*v(2), A1(2,2)*v(1)+A2(2,2)*v(2)]);


M_div_h2 = diag([m1,m2])/delta_t^2;
           
%% declare aux_vals
aux_vals = struct('M_div_h2', M_div_h2, ...
                  'g',  g, ...
                  'f_k', f_k, ...
                  'f_K', f_K, ...
                  'rhs', rhs, ...
                  'dK_dk1', dK_dk1, ...
                  'dK_dk2', dK_dk2, ...
                  'f_dk_dx', f_dk_dx, ...
                  'f_tens_vec', f_tens_vec, ...
                  'tol', tol);


%% equilibrium

[y2,k2,err] = newton_equi(x_init, aux_vals);
relErr = err([1 1:end-1]) ./ err;
order = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))))

figure;
plot(relErr);
title('Relative Error');
grid on;

disp('Equilibrium solution');
abs(-y2(1)*k2(1) - (m1+m2)*g) < tol
abs(-y2(2)*k2(2) - m2*g) < tol

%% simulation

[t,y,err1] = BW_2DOF_scheme(t, x_init, aux_vals);

%% plots
figure;
plot(t,y);
hold on;
y_equi = y2 * ones(1, length(y(1,:)));
plot(t,y_equi);
title('Position');
grid on;
legend('x1','x2','x1 equi', 'x2 equi');

figure;
plot(t,err1);
title('Absolute Error');


