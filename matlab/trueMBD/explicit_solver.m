function [t,x_vector_new] = explicit_solver(f, t, x_previous)
    
%     [~, x_vector_new] = ode45(f, [t(1), t(2)], x_vector);


    %initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];
    
    for n = 2 : length(t)
        
        delta_t = t(n) - t(n-1);
        x_previous = x_vector_new(n-1, :);
        f_n = f(t(n-1), x_previous')';
        x_vector_new(n,:) = x_previous + delta_t * f_n;

    end
    
      
end
