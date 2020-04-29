function [qc, pw1, pw2, pw3, pw4, pt1, pt2, pt3, pt4] = get_initial_positions(qc, r1, r2, r3, r4, pcc, upper_length, lower_length)
        qc = qc/norm(qc);

        C_Nc = get_basis(qc);      % get the transform from car to normal basis

        global_y = C_Nc(:,2);       % direction of the legs
        global_y = -global_y / norm(global_y);
        
        global_r1 = pcc + C_Nc*r1;        % global positions of the car corners
        global_r2 = pcc + C_Nc*r2;
        global_r3 = pcc + C_Nc*r3;
        global_r4 = pcc + C_Nc*r4;
        
        upper_global_spring_1 = upper_length(1)*global_y;   % elongation vectors of the spring
        upper_global_spring_2 = upper_length(2)*global_y;
        upper_global_spring_3 = upper_length(3)*global_y;
        upper_global_spring_4 = upper_length(4)*global_y;
        
        lower_global_spring_1 = lower_length(1)*global_y;  
        lower_global_spring_2 = lower_length(2)*global_y;
        lower_global_spring_3 = lower_length(3)*global_y;
        lower_global_spring_4 = lower_length(4)*global_y;
        
        pw1 = global_r1 + upper_global_spring_1;
        pw2 = global_r2 + upper_global_spring_2;
        pw3 = global_r3 + upper_global_spring_3;
        pw4 = global_r4 + upper_global_spring_4;

        pt1 = pw1 + lower_global_spring_1;
        pt2 = pw2 + lower_global_spring_2;
        pt3 = pw3 + lower_global_spring_3;
        pt4 = pw4 + lower_global_spring_4;
end