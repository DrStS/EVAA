function [x_curr,k] =  newton_equi(x_curr, aux_vals)
% get aux_vals
m1 = aux_vals.m1;
m2 = aux_vals.m2;
g = aux_vals.g;
l1 = aux_vals.l1;
l2 = aux_vals.l2;
tol = aux_vals.tol;

% define functions to update variables
f_k = aux_vals.interpolate_k;
f_rhs = @(k)([-m1*g + k(1)*l1 - k(2)*l2; ...
         -m2*g + k(2)*l2]);
f_K = @(k)([k(1)+k(2), -k(2); -k(2), k(2)]);

% fuction for newton loop
f_newton = @(x,K,rhs)(K*x - rhs);

% initial values for K and right hand side
k = f_k(x_curr);
K = f_K(k);
rhs = f_rhs(k);
    
while 1
    %dK/dk1
    dK1 = [1, 0; 0, 0];
    %dK/dk2
    dK2 = [1, -1; -1, 1];
    %dk1/dx
    dk1 = [2*(3*x_curr(1)- x_curr(2)- 2*l1+ l2), 2*(x_curr(2)-x_curr(1)-l2)];
    %dk2/dx
    dk2 = [2*(6*x_curr(1)- 5*x_curr(2)- l1 + 5*l2), 10*(x_curr(2)-x_curr(1)-l2)];
    %dk/dx
    dk = [dk1; dk2];
    %drhs/dk
    drhs = [l1, -l2; 0, l2];
    % J = K + dK/dk1 * (dk1/dx * x) + dK/dk2 * (dk2/dx * x) - drhs/dk * dk/dx 
    J = K + dK1*(dk1*x_curr)+dK2*(dk2*x_curr)-drhs*dk;
        
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