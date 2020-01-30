function [t,y] = BW_2DOF_scheme(t, x_prev, x_curr, aux_vals)

% get auxiliary variables
y = zeros(2, length(t));
y(:,1) = x_prev;
y(:,2) = x_curr;
h = t(2)-t(1);

m1 = aux_vals.m1;
m2 = aux_vals.m2;
g = aux_vals.g;
l1 = aux_vals.l1;
l2 = aux_vals.l2;
tol = aux_vals.tol;

% calc not changing matrices
M = diag([aux_vals.m1, aux_vals.m2]);
M_div_h2 = M / h^2;


% define functions to update variables
f_k = aux_vals.interpolate_k;
f_rhs = @(k)([-m1*g + k(1)*l1 - k(2)*l2; ...
         -m2*g + k(2)*l2]);
f_K = @(k)([k(1)+k(2), -k(2); -k(2), k(2)]);

% fuction for newton loop
f_newton = @(y,i,K,rhs)((M_div_h2 + K) * y(:, i+1)- 2*M_div_h2*y(:,i)+ M_div_h2*y(:,i-1)-rhs);

% init vars
k = f_k(y(:,2));
K = f_K(k);
rhs = f_rhs(k);

for i = 2 : length(t)-1
    
    % Create the linear system - init step
    y(:,i+1) = (M_div_h2 + K) \ ((2 * M_div_h2)*y(:,i) - M_div_h2*y(:,i-1)+rhs);
    
    %update values
    k = f_k(y(:,i+1));
    K = f_K(k);
    rhs = f_rhs(k);
    
    % newton loop
    while 1
        %dK/dk1
        dK1 = [1, 0; 0, 0];
        %dK/dk2
        dK2 = [1, -1; -1, 1];
        %dk1/dx
        dk1 = [2*(3*y(1,i+1)- y(2,i+1)- 2*l1+ l2), 2*(y(2,i+1)-y(1,i+1)-l2)];
        %dk2/dx
        dk2 = [2*(6*y(1,i+1)- 5*y(2,i+1)- l1 + 5*l2), 10*(y(2,i+1)-y(1,i+1)-l2)];
        %dk/dx
        dk = [dk1; dk2];
        %drhs/dk
        drhs = [l1, -l2; 0, l2];
        % J = M/h^2 + K + dK/dk1 * (dk1/dx * x) + dK/dk2 * (dk2/dx * x) - drhs/dk * dk/dx
        J = M_div_h2 + K + dK1*(dk1*y(:,i+1))+dK2*(dk2*y(:,i+1))-drhs*dk;
        
        %newton step
        y(:,i+1) = J\(-f_newton(y,i,K,rhs)+J*y(:,i+1));
        
        %update values to check error
        k = f_k(y(:,i+1));
        K = f_K(k);
        rhs = f_rhs(k);
        
        % if error < tol => break
        if (norm(f_newton(y,i,K,rhs)) < tol)
            break;
        end
    end  
end

end