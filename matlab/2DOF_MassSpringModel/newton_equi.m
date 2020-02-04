function [x_curr,k,err] =  newton_equi(x_curr, aux_vals)
% get aux_vals
tol = aux_vals.tol;

dK1 = aux_vals.dK_dk1; % matrix
dK2 = aux_vals.dK_dk2; % matrix
drhs = aux_vals.df_rhs_dk; % matrix

% define functions to update variables
f_k = aux_vals.interpolate_k; % function
f_K = aux_vals.f_K; % function
f_rhs = aux_vals.f_rhs; % function
f_dk_ddx = aux_vals.f_dk_ddx; % function
ddx_dx = aux_vals.ddx_dx;
get_dx = aux_vals.get_dx; %function

% fuction for newton loop
f_newton = @(x,K,rhs)(K*x - rhs);

% initial values for K and right hand side
k = f_k(get_dx(x_curr));
K = f_K(k);
rhs = f_rhs(k);

iter = 0;
err = [];
while 1
    iter = iter +1;
    
    dx = get_dx(x_curr);
    %dk/ddx
    dk_ddx = f_dk_ddx(dx);
    % J = K + dK/dk1*dk1/ddx * ddx/dx * x(1) + dK/dk2 * dk2/ddx * ddx/dx * x(2) - drhs/dk * dk/ddx * ddx/dx 
    J = K + dK1*dk_ddx*ddx_dx*x_curr(1)+dK2*dk_ddx*ddx_dx*x_curr(2)-drhs*dk_ddx*ddx_dx;
    
    %newton step
    dx = - J\f_newton(x_curr,K,rhs);
    x_curr = x_curr + dx;
    
    %update k, K and rhs to check the error;
    k = f_k(get_dx(x_curr));
    K = f_K(k);
    rhs = f_rhs(k);
    
    err = [err, norm(f_newton(x_curr,K,rhs),1)];
    % if error < tol => break
    if (norm(f_newton(x_curr,K,rhs),1) < tol)
        break;
    end
end
err = [err, norm(f_newton(x_curr,K,rhs))];
iter
end