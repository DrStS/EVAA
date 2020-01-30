function [t,y] = BW_2DOF_scheme(t, x_prev, x_curr, aux_vals)

%% Extract data
% get auxiliary variables
y = zeros(2, length(t));
y(:,1) = x_prev;
y(:,2) = x_curr;
h = t(2)-t(1);

m1 = aux_vals.m1;
m2 = aux_vals.m2;
g = aux_vals.g;
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

% calc not changing matrices
M = diag([m1, m2]);
M_div_h2 = M / h^2;

%% fuction for newton loop
f_newton = @(y,i,K,rhs)((M_div_h2 + K) * y(:, i+1)- 2*M_div_h2*y(:,i)+ M_div_h2*y(:,i-1)-rhs);

% init vars
k = f_k(get_dx(y(:,2)));
K = f_K(k);
rhs = f_rhs(k);

for i = 2 : length(t)-1
    
    % Create the linear system - init step
    y(:,i+1) = (M_div_h2 + K) \ ((2 * M_div_h2)*y(:,i) - M_div_h2*y(:,i-1)+rhs);
    
    %update values
    k = f_k(get_dx(y(:,i+1)));
    K = f_K(k);
    rhs = f_rhs(k);
    
    % newton loop
    while 1
        dx = get_dx(y(:,i+1));
        %dk/dx
        dk = dk_d_dx(dx);
        % J = M/h^2 + K + dK/dk1 * (dk1/dx * x) + dK/dk2 * (dk2/dx * x) - drhs/dk * dk/dx
        J = M_div_h2 + K + dK1*(dk(1,:)*ddx_dx*y(:,i+1))+dK2*(dk(2,:)*ddx_dx*y(:,i+1))-drhs*dk;
        
        %newton step
        y(:,i+1) = J\(-f_newton(y,i,K,rhs)+J*y(:,i+1));
        
        %update values to check error
        k = f_k(dx);
        K = f_K(k);
        rhs = f_rhs(k);
        
        % if error < tol => break
        if (norm(f_newton(y,i,K,rhs)) < tol)
            break;
        end
    end  
end

end