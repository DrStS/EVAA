function [] = error_plotter3D(  solver, t, num_iter, y, y_sol, mass, mass_wheel, mass_tyre, I, g,...
                                upper_spring_length, lower_spring_length, ...
                                upper_spring_stiffness, lower_spring_stiffness, ...
                                upper_rotational_stiffness, lower_rotational_stiffness)
    figure()
   
    %% error measure
    error = zeros(size(t,2), 1);
                      
    for i = 1 : num_iter

        %% get components
        pcc = y_sol(i,35:37)';

        qc = y_sol(i,31:34)';
        
        basis_c = get_basis(qc);
    
        % apply C_Cos transform
        C_cN = C_cos_transf(basis_c, eye(3));

        pc1 = y(i,8:10)';
        pc2 = y(i,11:13)';
        pc3 = y(i,14:16)';
        pc4 = y(i,17:19)';

        
        pw1 = y_sol(i,38:40)';
        pw2 = y_sol(i,41:43)';
        pw3 = y_sol(i,44:46)';
        pw4 = y_sol(i,47:49)';
        pw = [pw1, pw2, pw3, pw4];
        
        pt1 = y_sol(i,50:52)';
        pt2 = y_sol(i,53:55)';
        pt3 = y_sol(i,56:58)';
        pt4 = y_sol(i,59:61)';
        pt = [pt1, pt2, pt3, pt4];
        
        r_up1 = pc1 - pw1;
        r_up2 = pc2 - pw2;
        r_up3 = pc3 - pw3;
        r_up4 = pc4 - pw4;

        r_low1 = pw1 - pt1;
        r_low2 = pw2 - pt2;
        r_low3 = pw3 - pt3;
        r_low4 = pw4 - pt4;

        %% get angle and normal vectors at the legs

        [~, upper_angle1, ~] = get_quaternion(r_up1, C_cN(2,:)');
        [~, upper_angle2, ~] = get_quaternion(r_up2, C_cN(2,:)');
        [~, upper_angle3, ~] = get_quaternion(r_up3, C_cN(2,:)');
        [~, upper_angle4, ~] = get_quaternion(r_up4, C_cN(2,:)');
        upper_angle = [upper_angle1, upper_angle2, upper_angle3, upper_angle4];

        [~, lower_angle1, ~] = get_quaternion(r_low1, r_up1);
        [~, lower_angle2, ~] = get_quaternion(r_low2, r_up2);
        [~, lower_angle3, ~] = get_quaternion(r_low3, r_up3);
        [~, lower_angle4, ~] = get_quaternion(r_low4, r_up4);
        lower_angle = [lower_angle1, lower_angle2, lower_angle3, lower_angle4];

        
        upper_length1 = norm(r_up1);
        upper_length2 = norm(r_up2);
        upper_length3 = norm(r_up3);
        upper_length4 = norm(r_up4);
        upper_length = [upper_length1, upper_length2, upper_length3, upper_length4];
        
        lower_length1 = norm(r_low1);
        lower_length2 = norm(r_low2);
        lower_length3 = norm(r_low3);
        lower_length4 = norm(r_low4);
        lower_length = [lower_length1, lower_length2, lower_length3, lower_length4];

        wc = y_sol(i,1:3)';

        vc = y_sol(i,4:6)';

        vw1 = y_sol(i,7:9)';
        vw2 = y_sol(i,10:12)';
        vw3 = y_sol(i,13:15)';
        vw4 = y_sol(i,16:18)';
        vw = [vw1, vw2, vw3, vw4];

        vt1 = y_sol(i,19:21)';
        vt2 = y_sol(i,22:24)';
        vt3 = y_sol(i,25:27)';
        vt4 = y_sol(i,28:30)';
        vt = [vt1, vt2, vt3, vt4];

        %compute the energies in the system
        potential_energy = g * (mass * pcc(2) + sum(mass_wheel .* pw(2,:)) + sum(mass_tyre .* pt(2,:)));
        kinetic_energy = 0.5 * (mass * sum(vc.^2) + sum(mass_wheel .* sum(vw.^2)) + sum(mass_tyre .* sum(vt.^2)));
        rotational_energy = 0.5 * wc' * I * wc;
        spring_potential = 0.5 * (sum(upper_spring_stiffness .* (upper_spring_length - upper_length').^2) ...
                                + sum(lower_spring_stiffness .* (lower_spring_length - lower_length').^2));
        rot_spring_potential =  0.5 * sum(upper_rotational_stiffness' .* upper_angle.^2) + ...
                                0.5 * sum(lower_rotational_stiffness' .* lower_angle.^2);
                            
        %compute the total energy
%        error(i) = potential_energy + kinetic_energy + rotational_energy + spring_potential + rot_spring_potential;
        error(i) = rot_spring_potential;

        % compute the relative error
        if i > 1
%            error(i) = (error(i) - error(1)) / error(1);
            error(i) = (error(i) - error(1));

        end
    end
        error(1) = 0;

    %% plot characteristics    
    plot(t(2:end-1), error(2:end-1)); grid on;
    hold on;
        
    legend('energy variation in the 3D system');
    title_array = ['Error graph - BDF2, ', num2str(t(end)), ' sec, dt=',num2str(t(2)-t(1))];
    title(title_array);
    xlabel('Time [sec]');
    ylabel('Energy variation');
end