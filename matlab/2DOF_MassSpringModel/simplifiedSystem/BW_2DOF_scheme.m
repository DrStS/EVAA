function [t,y,err] = BW_2DOF_scheme(t, x_curr, aux_vals)

%% Extract data
% get auxiliary variables
M_div_h2 = aux_vals.M_div_h2;
g = aux_vals.g;
tol = aux_vals.tol;
dK1 = aux_vals.dK_dk1; % matrix
dK2 = aux_vals.dK_dk2; % matrix
rhs = aux_vals.rhs; % matrix

% tensor vector multiplication (2,2,2)x(2,1) -> (2,2)
f_tens_vec = aux_vals.f_tens_vec;

y = zeros(2, length(t));
err = zeros(length(t));
y(:,1) = x_curr;
y(:,2) = x_curr;

% define functions to update variables
f_k = aux_vals.f_k; % function
f_K = aux_vals.f_K; % function
f_dk_dx = aux_vals.f_dk_dx; % function

%% fuction for newton loop
f_newton = @(y,i,K)((M_div_h2 + K) * y(:, i+1)- 2*M_div_h2*y(:,i)+ M_div_h2*y(:,i-1)-rhs);

% init vars
k = f_k(x_curr);
K = f_K(k);

for i = 2 : length(t)-1
    % newton loop
    j = 0;
    while 1
        j = j + 1;
        %dk/dx
        dk_dx = f_dk_dx(y(:,i+1));
        % J = M/h^2 + K + (dK/dk1 * dk1/ddx, dK/dk2 * dk2/ddx) * x
        J = M_div_h2 + K + f_tens_vec(dK1*dk_dx,dK2*dk_dx,y(:,i+1));
        %newton step
        y(:,i+1) = J\(J*y(:,i+1)-f_newton(y,i,K));
        
        %update values to check error
        k = f_k(y(:,i+1));
        K = f_K(k);
        
        % if error < tol => break
        err(i) = norm(f_newton(y,i,K));
        if (err(i) < tol | j == 15)
            break;
        end
    end  
end

end