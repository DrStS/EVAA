function [t, x_vector_new, metrics] = Broyden_Euler(f, t, x_previous, tol, max_iter)
    % Initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];
     metrics = zeros(length(t), 2); %first component: number of iteration to convergence, second component: condition of the Jacobian
   
    dt_inv = round(1/(t(2)-t(1)));
    
    for n = 2 : length(t)
        delta_t = t(n) - t(n-1);
        % 0 ~= F(x_new) = x_new - x_n - delta_t * x_dot(x_new)
        
        x_previous = x_vector_new(n-1, :);
        
        % 1. Initialize guess from previous time step
        f_old = f(t(n-1),n-1, x_previous');
        % in case the velocity is 0 add nuggets to avoid singular matrices
        nuggets = 1e-4;
        f_old(abs(f_old) < nuggets) = nuggets*sign(f_old(abs(f_old) < nuggets));
        f_old(f_old==0) = (2*randi(2)-3)*nuggets;
       
        % 2. Initial guess from Explicit Euler
        x = x_previous + delta_t * f_old';
        f_new = f(t(n),n, x');
        
        % Initial approximation of the Jacobian
        dx = x - x_previous;
        df = f_new' - f_old';
        
        % approximate J(x_0)         
        J = eye(length(df)) - delta_t*((1./dx)'*df)'; 
        metrics(n,2) = cond(J);
        
        % calculate initial F for stopping condition
        x_dot = f_new';
        F = dx - delta_t * x_dot;
        
        % Broyden's Method
        for i = 1 : max_iter
            if (norm(F) < tol)
                break;
            end
            
            % x(i+1) = x(i) - J^(-1)*g(x(i))
            x_new = (x' - J\F')';
            
            % Calculate new derivative
            x_dot = f(t(n-1),n, x_new')';
            
            F_new = x_new - x_previous - delta_t * x_dot;
            
            % J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
            dF = (F_new - F)';
            dx = (x_new - x)';
            J = J + ((dF/ norm(dx)) - J * (dx/ norm(dx))) * (dx / norm(dx))'; 

            F = F_new;
            x = x_new;
        end
        metrics(n,1) = i;
        if(i==max_iter)
            dispstr = ['Maximum number of iteration ', num2str(max_iter),' reached! Current accuracy: ', num2str(norm(F))];
            disp(dispstr)
        end
        if (mod(n, dt_inv)==0)
            timestr = ['Iteration ', num2str(n), ' at time ', num2str(t(n+1))];
            disp(timestr);
        end
		
        x_vector_new(n,:) = x;
    end
end