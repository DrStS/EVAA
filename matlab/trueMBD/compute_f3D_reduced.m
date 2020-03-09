function f = compute_f3D_reduced(x_vector,t,aux_vals)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%                 Setup                  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% extract data from x_vector
    wc = x_vector(1:3);
    vc = x_vector(4:6);
    vw1 = x_vector(7:9);
    vw2 = x_vector(10:12);
    vw3 = x_vector(13:15);
    vw4 = x_vector(16:18);
    vt1 = x_vector(19:21);
    vt2 = x_vector(22:24);
    vt3 = x_vector(25:27);
    vt4 = x_vector(28:30);
    qc = x_vector(31:34);
    pcc = x_vector(35:37);
    pw1 = x_vector(38:40);
    pw2 = x_vector(41:43);
    pw3 = x_vector(44:46);
    pw4 = x_vector(47:49);
    pt1 = x_vector(50:52);
    pt2 = x_vector(53:55);
    pt3 = x_vector(56:58);
    pt4 = x_vector(59:61);
    

    %% extract data from aux_vals
    r1 = aux_vals.r1;                                           %in local car basis
    r2 = aux_vals.r2;
    r3 = aux_vals.r3;
    r4 = aux_vals.r4;
    r1_tilda = aux_vals.r1_tilda;                               %in local car basis
    r2_tilda = aux_vals.r2_tilda;
    r3_tilda = aux_vals.r3_tilda;
    r4_tilda = aux_vals.r4_tilda;
    FC = [0; aux_vals.FC; 0];                                   %in global basis
    FW = aux_vals.FW;                                           %in global basis
    FT = aux_vals.FT;                                           %in global basis    
    lower_spring_length = aux_vals.lower_spring_length;         %in local leg basis
    upper_spring_length = aux_vals.upper_spring_length;
    lower_spring_stiffness = aux_vals.lower_spring_stiffness;
    upper_spring_stiffness = aux_vals.upper_spring_stiffness;
    lower_rotational_stiffness = aux_vals.lower_rotational_stiffness;
    upper_rotational_stiffness = aux_vals.upper_rotational_stiffness;
    
    A = aux_vals.A;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%                 Basis                  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% get cosine transforms (C_Nc means r_N = C_Nc * r_c)
    % compute local base vectors
    basis_N = eye(3);

    basis_c = get_basis(qc);
    
    % apply C_Cos transform
    C_cN = C_cos_transf(basis_c, basis_N);      %no need to apply it, should be C_cN = basis_c'

    % global positions of the tyre connections
    pc1 = C_cN' * r1 + pcc;
    pc2 = C_cN' * r2 + pcc;
    pc3 = C_cN' * r3 + pcc;
    pc4 = C_cN' * r4 + pcc;

    % currrent length of the legs
    r_up1 = pc1 - pw1;
    r_up2 = pc2 - pw2;
    r_up3 = pc3 - pw3;
    r_up4 = pc4 - pw4;

    r_low1 = pw1 - pt1;
    r_low2 = pw2 - pt2;
    r_low3 = pw3 - pt3;
    r_low4 = pw4 - pt4;

    %% get angle and normal vectors at the legs

    [~, upper_angle1, upper_normal1] = get_quaternion(r_up1, C_cN(2,:)');
    [~, upper_angle2, upper_normal2] = get_quaternion(r_up2, C_cN(2,:)');
    [~, upper_angle3, upper_normal3] = get_quaternion(r_up3, C_cN(2,:)');
    [~, upper_angle4, upper_normal4] = get_quaternion(r_up4, C_cN(2,:)');

    [~, lower_angle1, lower_normal1] = get_quaternion(r_low1, r_up1);
    [~, lower_angle2, lower_normal2] = get_quaternion(r_low2, r_up2);
    [~, lower_angle3, lower_normal3] = get_quaternion(r_low3, r_up3);
    [~, lower_angle4, lower_normal4] = get_quaternion(r_low4, r_up4);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%           Forces  and Torques          %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% calculate the elongational spring forces (in global basis)
    upper_force1 = upper_spring_stiffness(1) * (r_up1) * (1 - upper_spring_length(1) / norm(r_up1));    
    upper_force2 = upper_spring_stiffness(2) * (r_up2) * (1 - upper_spring_length(2) / norm(r_up2));
    upper_force3 = upper_spring_stiffness(3) * (r_up3) * (1 - upper_spring_length(3) / norm(r_up3));
    upper_force4 = upper_spring_stiffness(4) * (r_up4) * (1 - upper_spring_length(4) / norm(r_up4));

    lower_force1 = lower_spring_stiffness(1) * (r_low1) * (1 - lower_spring_length(1) / norm(r_low1));
    lower_force2 = lower_spring_stiffness(2) * (r_low2) * (1 - lower_spring_length(2) / norm(r_low2));
    lower_force3 = lower_spring_stiffness(3) * (r_low3) * (1 - lower_spring_length(3) / norm(r_low3));
    lower_force4 = lower_spring_stiffness(4) * (r_low4) * (1 - lower_spring_length(4) / norm(r_low4));

    
    %% torque from the rotational spring
    upper_S1 = upper_rotational_stiffness(1) * upper_angle1 * upper_normal1;        % in global basis
    upper_S2 = upper_rotational_stiffness(2) * upper_angle2 * upper_normal2;
    upper_S3 = upper_rotational_stiffness(3) * upper_angle3 * upper_normal3;
    upper_S4 = upper_rotational_stiffness(4) * upper_angle4 * upper_normal4; 

    lower_S1 = lower_rotational_stiffness(1) * lower_angle1 * lower_normal1;        % in global basis
    lower_S2 = lower_rotational_stiffness(2) * lower_angle2 * lower_normal2;
    lower_S3 = lower_rotational_stiffness(3) * lower_angle3 * lower_normal3;
    lower_S4 = lower_rotational_stiffness(4) * lower_angle4 * lower_normal4; 


    %% calculate the effect of the rotational springs (in global basis)
    lower_rot_force1 = -cross( lower_S1, r_low1) / (r_low1'*r_low1);     
    lower_rot_force2 = -cross( lower_S2, r_low2) / (r_low2'*r_low2);
    lower_rot_force3 = -cross( lower_S3, r_low3) / (r_low3'*r_low3);
    lower_rot_force4 = -cross( lower_S4, r_low4) / (r_low4'*r_low4);

    upper_rot_force1 = -cross( upper_S1, r_up1) / (r_up1'*r_up1);
    upper_rot_force2 = -cross( upper_S2, r_up2) / (r_up2'*r_up2);
    upper_rot_force3 = -cross( upper_S3, r_up3) / (r_up3'*r_up3);
    upper_rot_force4 = -cross( upper_S4, r_up4) / (r_up4'*r_up4);

    car_rot_force1 = -cross( lower_S1, r_up1) / (r_up1'*r_up1);
    car_rot_force2 = -cross( lower_S2, r_up2) / (r_up2'*r_up2);
    car_rot_force3 = -cross( lower_S3, r_up3) / (r_up3'*r_up3);
    car_rot_force4 = -cross( lower_S4, r_up4) / (r_up4'*r_up4);
    
    car_rot_force = car_rot_force1 + car_rot_force2 + car_rot_force3 + car_rot_force4;
    
    %% get external forces    
    % for the legs
    local_FW1 = [0; FW(1); 0];      %in global basis
    local_FW2 = [0; FW(2); 0];
    local_FW3 = [0; FW(3); 0];
    local_FW4 = [0; FW(4); 0];
    
    local_FT1 = [0; FT(1); 0];      %in global basis
    local_FT2 = [0; FT(2); 0];
    local_FT3 = [0; FT(3); 0];
    local_FT4 = [0; FT(4); 0];

    % road forces    
    local_FR1 = [0;  aux_vals.FR1(t, pt1(2), vt1(2), A(19,19), lower_force1(2) + FT(1) + lower_rot_force1(2)); 0];      %in global basis
    local_FR2 = [0;  aux_vals.FR2(t, pt2(2), vt2(2), A(22,22), lower_force2(2) + FT(2) + lower_rot_force2(2)); 0];   
    local_FR3 = [0;  aux_vals.FR3(t, pt3(2), vt3(2), A(25,25), lower_force3(2) + FT(3) + lower_rot_force3(2)); 0]; 
    local_FR4 = [0;  aux_vals.FR4(t, pt4(2), vt4(2), A(28,28), lower_force4(2) + FT(4) + lower_rot_force4(2)); 0]; 


    %% get H=I*w
    Hc = A(1:3, 1:3) * wc;           %in local car basis

    %% external torque on the car body (use later for rotational damping) 
    Tc = zeros(3,1);     %in local car basis

    
    %% sum of all torques induced by forces
    % for the car body
    sum_torque_spring_car = r1_tilda * (C_cN * (-upper_force1 + car_rot_force1)) + ...     % from the elongational springs
                            r2_tilda * (C_cN * (-upper_force2 + car_rot_force2)) + ...
                            r3_tilda * (C_cN * (-upper_force3 + car_rot_force3)) + ...
                            r4_tilda * (C_cN * (-upper_force4 + car_rot_force4)) + ...
                           -C_cN * upper_S1 - C_cN * upper_S2 - C_cN * upper_S3 - C_cN * upper_S4 + ...     % ??from the rotational spring
                            -get_tilda(wc) * Hc + Tc;                       % from angular momentum and external torques

                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%                 Solve                  %%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    b =  [sum_torque_spring_car;...                                                           %w_dot_c
         FC - upper_force1 - upper_force2 - upper_force3 - upper_force4 + car_rot_force + ...
              - upper_rot_force1 - upper_rot_force2 - upper_rot_force3 - upper_rot_force4; ...%vc_dot                                                         
         upper_force1 - lower_force1 + local_FW1 + upper_rot_force1 + ...
              - car_rot_force1 - lower_rot_force1; ...                                        %vw1_dot
         upper_force2 - lower_force2 + local_FW2 + upper_rot_force2 + ...
              - car_rot_force2 - lower_rot_force2; ...                                        %vw2_dot
         upper_force3 - lower_force3 + local_FW3 + upper_rot_force3 + ...
              - car_rot_force3 - lower_rot_force3; ...                                        %vw3_dot
         upper_force4 - lower_force4 + local_FW4 + upper_rot_force4 + ...
              - car_rot_force4 - lower_rot_force4; ...                                        %vw4_dot
         lower_force1 + local_FT1 + local_FR1 + lower_rot_force1; ...                         %vt1_dot
         lower_force2 + local_FT2 + local_FR2 + lower_rot_force2; ...                         %vt2_dot
         lower_force3 + local_FT3 + local_FR3 + lower_rot_force3; ...                         %vt3_dot
         lower_force4 + local_FT4 + local_FR4 + lower_rot_force4];                            %vt4_dot
     
     % solve the system
    result_vector = A \ b;
    
    
    %% get the derivative of the attitude (expressed in quaternions) from the angular velocities
    Qc = 0.5 * [qc(4) -qc(3) qc(2); qc(3) qc(4) -qc(1); -qc(2) qc(1) qc(4); -qc(1) -qc(2) -qc(3)];
    qc_dot = Qc * wc;

    
    %% evaluate the function -> x_dot = f(x)
    f = [result_vector(1:3); ...    % wc_dot
         result_vector(4:6); ...    % vc_dot
         result_vector(7:9); ...    % vw1_dot
         result_vector(10:12); ...  % vw2_dot
         result_vector(13:15); ...  % vw3_dot
         result_vector(16:18); ...  % vw4_dot
         result_vector(19:21); ...  % vt1_dot
         result_vector(22:24); ...  % vt2_dot
         result_vector(25:27); ...  % vt3_dot
         result_vector(28:30); ...  % vt4_dot
         qc_dot; ...                % qc_dot
         vc; ...                    % pcc_dot
         vw1; ...                   % pw1_dot
         vw2; ...                   % pw2_dot
         vw3; ...                   % pw3_dot
         vw4; ...                   % pw4_dot
         vt1; ...                   % pt1_dot
         vt2; ...                   % pt2_dot
         vt3; ...                   % pt3_dot
         vt4];                      % pt4_dot

end