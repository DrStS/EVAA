function [t,x_vector_new] = Runge_Kutta_4(f, t, x_previous)

    %initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];
    
    for n = 2 : length(t)

        delta_t = t(n) - t(n-1);
        t_prev = t(n-1);
        x_previous = x_vector_new(n-1, :);
        k1 = delta_t * f(t_prev, x_previous')';
        k2 = delta_t * f(t_prev + delta_t/2, x_previous' + k1'/2)';
        k3 = delta_t * f(t_prev + delta_t/2, x_previous' + k2'/2)';
        k4 = delta_t * f(t_prev + delta_t, x_previous' + k3')';
        
        x_vector_new(n,:) = x_previous + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

    end
    
      
end
