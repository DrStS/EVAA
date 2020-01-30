function [x_curr,k] =  newton_equi(x_curr, aux_vals)
% get aux_vals
tol = aux_vals.tol;

dK1 = aux_vals.dK_dk1; % matrix
dK2 = aux_vals.dK_dk2; % matrix
drhs = aux_vals.df_rhs_dk; % matrix

% define functions to update variables
f_k = aux_vals.interpolate_k; % function
f_K = aux_vals.f_K; % function
f_rhs = aux_vals.f_rhs; % function
dk_d_dx = aux_vals.dk_d_dx; % function
ddx_dx = aux_vals.ddx_dx; % function
get_dx = aux_vals.get_dx; % function

% fuction for newton loop
f_newton = @(x,K,rhs)(K*x - rhs);

% initial values for K and right hand side
k = f_k(x_curr);
K = f_K(k);
rhs = f_rhs(k);
    
while 1
    dx = get_dx(x_curr);
    %dk/dx
    dk = dk_d_dx(dx);
    % J = K + dK/dk1 * (dk1/dx * x) + dK/dk2 * (dk2/dx * x) - drhs/dk * dk/dx 
    J = K + dK1*(dk(1,:)*ddx_dx*x_curr)+dK2*(dk(2,:)*ddx_dx*x_curr)-drhs*dk;
        
    %newton step
    x_curr = J\(-f_newton(x_curr,K,rhs)+J*x_curr);
    
    %update k, K and rhs to check the error;
    k = f_k(x_curr);
    K = f_K(k);
    rhs = f_rhs(k);
    
    % if error < tol => break
    if (norm(f_newton(x_curr,K,rhs)) < tol)
        break;
    end
end
end