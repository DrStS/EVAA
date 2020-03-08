function [t, x_vector_new] = Broyden_Crank_Nicolson(f, t, x_previous, tol, max_iter)
    
    % Initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];
    
    for n = 2 : length(t)
		delta_t = t(n) - t(n-1);
		% 0 ~= F(x_new) = x_new - x_n - delta_t * x_dot(x_new)

		x_previous = x_vector_new(n-1, :);
		
		% 1. Initial guess using previous time step
		f_old = f(t(n-1), x_previous');
		% in case the velocity is 0 add nuggets to avoid singular matrices
		f_old(abs(f_old) < 0.001) = 0.001;
	   
		% 2. Initial guess from Explicit Euler
		x = x_previous + delta_t * f_old';
		f_new = f(t(n), x');
		
		% Initial approximation of the Jacobian
		dx = x - x_previous;
		df = f_new' - f_old';
		
		% Approximate J(x_0)         
		J = eye(length(df)) - delta_t *0.5 * ((1./dx)' * df)'; 
		
		% Calculate initial F for stopping condition
		x_dot = f_new';
		x_dot_previous = f_old';
		F = dx - delta_t * 0.5 * (x_dot + x_dot_previous);
		
		% Broyden's Method
		for i = 1 : max_iter
			if (norm(F) < tol)
				break;
			end
			
			% x(i+1) = x(i) - J^(-1)*g(x(i))
			x_new = (x' - J\F')';
			
			% Calculate new derivative
			x_dot = f(t(n-1), x_new')';
			
			F_new = x_new - x_previous - delta_t * 0.5 * (x_dot + x_dot_previous);
			
			% J(n+1) = J(n) + (dF - J * dx) * dx' / norm(dx)^2
			dF = (F_new - F)';
			dx = (x_new - x)';
			J = J + ((dF/ norm(dx)) - J * (dx/ norm(dx))) * (dx' / norm(dx)); 

			F = F_new;
			x = x_new;
		end
		x_vector_new(n,:) = x;
    end  
end