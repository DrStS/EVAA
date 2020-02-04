function [x_curr,k,err] =  newton_equi(x_curr, aux_vals)
% get aux_vals
tol = aux_vals.tol;

rhs = aux_vals.rhs;
dK1 = aux_vals.dK_dk1; % matrix
dK2 = aux_vals.dK_dk2; % matrix

% define functions to update variables
f_k = aux_vals.f_k; % function
f_K = aux_vals.f_K; % function
f_dk_dx = aux_vals.f_dk_dx; % function

% tensor vector multiplication (2,2,2)x(2,1) -> (2,2)
f_tens_vec = aux_vals.f_tens_vec;

% fuction for newton loop
f_newton = @(x,K)(K*x - rhs);

% initial values for K and right hand side
k = f_k(x_curr);
K = f_K(k);

iter = 0;
err = [];
while 1
    iter = iter +1;
    disp('iter');
    %dk/ddx
    dk_dx = f_dk_dx(x_curr);
    % J = K + (dK/dk1 * dk1/ddx, dK/dk2 * dk2/ddx) * x
    J = K + f_tens_vec(dK1*dk_dx,dK2*dk_dx,x_curr);

    
    err = [err, norm(f_newton(x_curr,K))];
    %newton step
    x_curr = J\(J*x_curr-f_newton(x_curr,K));
    
    %update k, K and rhs to check the error;
    k = f_k(x_curr);
    K = f_K(k);
    
    % if error < tol => break
    if (norm(f_newton(x_curr,K)) < tol)
        break;
    end
end
err = [err, norm(f_newton(x_curr,K))];
iter
end