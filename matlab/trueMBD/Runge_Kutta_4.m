function [t,x_vector_new, metrics] = Runge_Kutta_4(f, t, x_previous)
    metrics = 0;
    
    %initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];
    
    dt_inv = round(1/(t(2)-t(1))); %only needed to print progress

    for n = 2 : length(t)

        delta_t = t(n) - t(n-1);
        t_prev = t(n-1);
        x_previous = x_vector_new(n-1, :);
        k1 = delta_t * f(t_prev, n, x_previous')';
        k2 = delta_t * f(t_prev + delta_t/2, n, x_previous' + k1'/2)';
        k3 = delta_t * f(t_prev + delta_t/2, n, x_previous' + k2'/2)';
        k4 = delta_t * f(t_prev + delta_t, n, x_previous' + k3')';
        
        x_vector_new(n,:) = x_previous + 1/6 * (k1 + 2*k2 + 2*k3 + k4);

        if (mod(n, dt_inv)==0)
            timestr = ['Iteration ', num2str(n), ' at time ', num2str(t(n+1))];
            disp(timestr);
        end

    end
    
      
end
