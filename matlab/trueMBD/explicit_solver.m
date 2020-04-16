function [t,x_vector_new, metrics] = explicit_solver(f, t, x_previous)

    global previous_solution_vector previous_force_vector internal_tyre_forces

    previous_force_vector = zeros(3,4);
    internal_tyre_forces = zeros(3,4);

    metrics =0;

    %initialize return vector
    x_vector_new = [x_previous'; zeros(length(t)-1, length(x_previous))];

    dt_inv = round(1/(t(2)-t(1))); %only needed to print progress    
    
    for n = 2 : length(t)
        delta_t = t(n) - t(n-1);
        x_previous = x_vector_new(n-1, :);

        previous_solution_vector = x_previous;
        
        f_n = f(t(n-1), n-1, x_previous')';

        previous_force_vector = internal_tyre_forces;
        
        x_vector_new(n,:) = x_previous + delta_t * f_n;
        if (mod(n, dt_inv)==0)
            timestr = ['Iteration ', num2str(n), ' at time ', num2str(t(n+1))];
            disp(timestr);
        end

    end
    
      
end
