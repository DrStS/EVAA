function [t,y,err] = BW_2DOF_scheme(t, x_curr, aux_vals)

%% Extract data
% get auxiliary variables
M_div_h2 = aux_vals.M_div_h2;
g = aux_vals.g;
tol = aux_vals.tol;
dK1 = aux_vals.dK_dk1; % matrix
dK2 = aux_vals.dK_dk2; % matrix
rhs = aux_vals.rhs; % matrix

y = zeros(2, length(t));
err = zeros(length(t));
y(:,1) = x_curr;
y(:,2) = x_curr;

% define functions to update variables
f_k = aux_vals.f_k; % function
f_K = aux_vals.f_K; % function
f_dk_dx = aux_vals.f_dk_dx; % function

%% fuction for newton loop
f_newton = @(y_curr,y1,y2,K)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);

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
        J = M_div_h2 + K + dK1*dk_dx * y(1,i+1) + dK2*dk_dx * y(2,i+1);
        %newton step
        dy = - J \ f_newton(y(:,i+1), y(:,i), y(:,i-1),K);
        y(:,i+1) = y(:,i+1)+ dy;
        
        %update values to check error
        k = f_k(y(:,i+1));
        K = f_K(k);
        
        % if error < tol => break
        error = norm(f_newton(y(:,i+1), y(:,i), y(:,i-1),K),1);
        err(i) = error;
        if (err(i) < tol || j == 10)
            break;
        end
    end  
end

end