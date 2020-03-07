function f = compute_f3D_reduced(x_vector,t,aux_vals)
    % extract data from x_vector
    wc = x_vector(1:3);
    vc = x_vector(4);
    vw = x_vector(5:8);
    vt = x_vector(9:12);
    qc = x_vector(13:16);
    pcc = x_vector(17);
    
    % extract data from aux_vals
    r1 = aux_vals.r1;
    r2 = aux_vals.r2;
    r3 = aux_vals.r3;
    r4 = aux_vals.r4;
    r1_tilda = aux_vals.r1_tilda;
    r2_tilda = aux_vals.r2_tilda;
    r3_tilda = aux_vals.r3_tilda;
    r4_tilda = aux_vals.r4_tilda;
    FC = aux_vals.FC;
    FW = aux_vals.FW;
    FT = aux_vals.FT;
    lower_spring_length = aux_vals.lower_spring_length;
    upper_spring_length = aux_vals.upper_spring_length;
    lower_spring_stiffness = aux_vals.lower_spring_stiffness;
    upper_spring_stiffness = aux_vals.upper_spring_stiffness;
    reduced = aux_vals.reduced;
    A = aux_vals.A;
    Ic = A(1:3, 1:3);

    
    %get local basis vectors using unit quaternion rotation 
    s = 1 / norm(qc)^2;

    basis_c = zeros(3);
    basis_c(:,1) = [1 - 2 * s * (qc(2)^2 + qc(3)^2); ...
                    2 * s * (qc(1) * qc(2) + qc(3) * qc(4)); ...
                    2 * s * (qc(1) * qc(3) - qc(2) * qc(4))];
                    
    basis_c(:,2) = [2 * s * (qc(1) * qc(2) - qc(3) * qc(4)); ...
                    1 - 2 * s * (qc(1)^2 + qc(3)^2); ...
                    2 * s * (qc(2) * qc(3) + qc(1) * qc(4))];
        
    basis_c(:,3) = [2 * s * (qc(1) * qc(3) + qc(2) * qc(4)); ...
                    2 * s * (qc(2) * qc(3) - qc(1) * qc(4)); ...
                    1 - 2 * s * (qc(1)^2 + qc(2)^2)];
                
    basis_N = eye(3);
        
    
    %get cosine transforms (C_Nc means r_N = C_Nc * r_c)
    C_Nc = C_cos_transf(basis_N, basis_c);    
      
    % get positions of the corners of the car body
    r1_global = C_Nc * r1;
    r2_global = C_Nc * r2;
    r3_global = C_Nc * r3;
    r4_global = C_Nc * r4;
    
    car_corners = pcc + [r1_global(2); r2_global(2); r3_global(2); r4_global(2)];
    
    % calculate the length of the springs
    upper_length = car_corners - x_vector(18:21);
    lower_length = x_vector(18:21) - x_vector(22:25);

    % calculate the spring forces (take opposite if top node is considered)
    upper_force = upper_spring_stiffness .* (upper_length - upper_spring_length);
    lower_force = lower_spring_stiffness .* (lower_length - lower_spring_length);
    
    if reduced
        lower_force = [0; 0; 0; 0];
    end
    
    % convert spring forces to local basis
    upper_F1 = C_Nc' * [0; -upper_force(1); 0];
    upper_F2 = C_Nc' * [0; -upper_force(2); 0];
    upper_F3 = C_Nc' * [0; -upper_force(3); 0];
    upper_F4 = C_Nc' * [0; -upper_force(4); 0];
    
    %get H
    Hc = Ic * wc;
   
    % external torque on the car body (use later for rotational damping)
    Tc = zeros(4,1);    
    
    % sum of all torques from the springs
    sum_torque = sum([r1_tilda * upper_F1, r2_tilda * upper_F2, r3_tilda * upper_F3, r4_tilda * upper_F4],2);

    % get road forces
    FR = [aux_vals.FR1(t,lower_force(1) + FT(1)); ...
        aux_vals.FR2(t,lower_force(2) + FT(2)); ...
        aux_vals.FR3(t,lower_force(3) + FT(3)); ...
        aux_vals.FR4(t,lower_force(4) + FT(4))];


    % get the right hand side
    b = [-get_tilda(wc) * Hc + sum(Tc) + sum_torque;...
         FC - sum(upper_force); ...
         upper_force - lower_force + FW; ...
         lower_force + FR + FT];
     
    b(2) = 0;       % prohibit rotation along the y-axis

    % solve the system
    result_vector = A \ b;
        
    % get the derivative of the attitude (expressed in quaternions) from
    % the angular velocities
    Qc = 0.5 * [qc(4) -qc(3) qc(2); qc(3) qc(4) -qc(1); -qc(2) qc(1) qc(4); -qc(1) -qc(2) -qc(3)];
    
    qc_dot = Qc * wc;
    

    % evaluate the function - x_dot = f(x)
    f = [result_vector(1:3); ...    % wc_dot
         result_vector(4); ...      % vc_dot
         result_vector(5:8); ...    % vw_dot
         result_vector(9:12); ...   % vt_dot
         qc_dot; ...                % qc_dot
         vc; ...                    % pcc_dot
         vw; ...                    % ps_dot
         vt];                       % pw_dot
     
end