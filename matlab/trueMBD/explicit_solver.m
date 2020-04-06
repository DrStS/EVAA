function [t,x_vector_new, metrics] = explicit_solver(f, t, x_previous)
    metrics =0;
%     [~, x_vector_new] = ode45(f, [t(1), t(2)], x_vector);


    %initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];

    dt_inv = round(1/(t(2)-t(1))); %only needed to print progress    
    
    for n = 2 : length(t)
        
        delta_t = t(n) - t(n-1);
        x_previous = x_vector_new(n-1, :);
        f_n = f(t(n-1), n-1, x_previous')';
        x_vector_new(n,:) = x_previous + delta_t * f_n;
        if (mod(n, dt_inv)==0)
            timestr = ['Iteration ', num2str(n), ' at time ', num2str(t(n+1))];
            disp(timestr);
        end

    end
    
      
end
