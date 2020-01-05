function [t, y_return, y] = main_nasa_car(r1, r2, r3, r4, mass, mass_wheel, mass_tyre, Ic, qc, reduced, ... 
        lower_spring_length, upper_spring_length, initial_lower_spring_length, initial_upper_spring_length, lower_spring_stiffness, upper_spring_stiffness,...
        vc, vw, vt, wc, FC, FW1, FW2, FW3, FW4, FT1, FT2, FT3, FT4, FR1, FR2, FR3, FR4,...
        num_iter, delta_t, solver)
 
    qc = qc / norm(qc);                %normalize the orientation
        
    pcc = [0; 0; 0];                    % initial position of the center of the car at the origin
    
    % initial local basis of the car body
    basis_c = zeros(3);
    basis_c(:,1) = [1 - 2 *(qc(2)^2 + qc(3)^2); ...
                    2 * (qc(1) * qc(2) + qc(3) * qc(4)); ...
                    2 * (qc(1) * qc(3) - qc(2) * qc(4))];

    basis_c(:,2) = [2 * (qc(1) * qc(2) - qc(3) * qc(4)); ...
                    1 - 2 * (qc(1)^2 + qc(3)^2); ...
                    2 * (qc(2) * qc(3) + qc(1) * qc(4))];

    basis_c(:,3) = [2 * (qc(1) * qc(3) + qc(2) * qc(4)); ...
                    2 * (qc(2) * qc(3) - qc(1) * qc(4)); ...
                    1 - 2 * (qc(1)^2 + qc(2)^2)];

    % change of basis matrices
    C_Nc = C_cos_transf(eye(3), basis_c);
    
    % positions of the lower corners of the car body
    pc1 = pcc + C_Nc * r1;
    pc2 = pcc + C_Nc * r2;
    pc3 = pcc + C_Nc * r3;
    pc4 = pcc + C_Nc * r4;
    
    % positions of the wheels
    pw1 = [pc1(1);  pc1(2)-initial_upper_spring_length(1); pc1(3)  ];
    pw2 = [pc2(1);  pc2(2)-initial_upper_spring_length(2); pc2(3)  ];
    pw3 = [pc3(1);  pc3(2)-initial_upper_spring_length(3); pc3(3)  ];
    pw4 = [pc4(1);  pc4(2)-initial_upper_spring_length(4); pc4(3)  ];
    
    % position of the tyres
    pt1 = [pc1(1);  pc1(2)-initial_upper_spring_length(1)-initial_lower_spring_length(1)    ; pc1(3)  ];
    pt2 = [pc2(1);  pc2(2)-initial_upper_spring_length(2)-initial_lower_spring_length(2)    ; pc2(3)  ];
    pt3 = [pc3(1);  pc3(2)-initial_upper_spring_length(3)-initial_lower_spring_length(3)    ; pc3(3)  ];
    pt4 = [pc4(1);  pc4(2)-initial_upper_spring_length(4)-initial_lower_spring_length(4)    ; pc4(3)  ];
    
    % computation of ri_tilda, ro_tilda
    r1_tilda = get_tilda(r1);
    r2_tilda = get_tilda(r2);
    r3_tilda = get_tilda(r3);
    r4_tilda = get_tilda(r4);

    % initialise solution vector
    x_vector = [wc; ...     % 3 
                vc; ...     % 1
                vw; ...     % 4
                vt; ...     % 4
                qc; ...     % 4
                pcc(2); ... % 1
                pw1(2); ... % 1
                pw2(2); ... % 1
                pw3(2); ... % 1
                pw4(2); ... % 1
                pt1(2); ... % 1
                pt2(2); ... % 1
                pt3(2); ... % 1
                pt4(2);];   % 1

            
    % create the reduced system matrix
    A = diag([1,1,1, mass, mass_wheel, mass_tyre]);  
    A(1:3, 1:3) = Ic;
    A(2, 1:3) = [0, 1, 0];          %prohibit rotation along the y-axis

    % create the structure with parameters used in functions' computation
    aux_vals = struct('r1', r1, ...
                      'r2', r2, ...
                      'r3', r3, ...
                      'r4', r4, ...
                      'r1_tilda', r1_tilda, ...
                      'r2_tilda', r2_tilda, ...
                      'r3_tilda', r3_tilda, ...
                      'r4_tilda', r4_tilda, ...
                      'FC', FC, ...
                      'FW', [FW1; FW2; FW3; FW4], ...
                      'FT', [FT1; FT2; FT3; FT4], ...
                      'FR1', FR1, ...
                      'FR2', FR2, ...
                      'FR3', FR3, ...
                      'FR4', FR4, ...
                      'lower_spring_length', lower_spring_length, ...
                      'upper_spring_length', upper_spring_length, ...
                      'lower_spring_stiffness', lower_spring_stiffness, ...
                      'upper_spring_stiffness', upper_spring_stiffness, ...
                      'reduced', reduced, ...
                      'A', A);


    t = 0:delta_t:num_iter*delta_t;

    f = @(t,x) compute_f3D_reduced(x,t,aux_vals);

    [t,y] = solver(f, t, x_vector);
   
    y_return = zeros(size(y,1), 13);  %components: qc(1:4), pcc(5), pc(6:9), pw(10:13), pt(14:17), x_coord(18:21), z_coord(22:25)

    y_return(:,1:5) = y(:,13:17);
    y_return(:,10:17) = y(:,18:25);

    % retrieve positions from the angles
    for i = 1:size(y,1)
        qc = y_return(i,1:4);
        pcc = [0; y_return(i,5); 0]; 
        
        %get local basis vectors using unit quaternion rotation 
        s = 1 / norm(qc)^2;

        basis_c = zeros(3);
        basis_c(:,1) = [1 - 2 * s *(qc(2)^2 + qc(3)^2); ...
                        2 * s * (qc(1) * qc(2) + qc(3) * qc(4)); ...
                        2 * s * (qc(1) * qc(3) - qc(2) * qc(4))];

        basis_c(:,2) = [2 * s * (qc(1) * qc(2) - qc(3) * qc(4)); ...
                        1 - 2 * s * (qc(1)^2 + qc(3)^2); ...
                        2 * s * (qc(2) * qc(3) + qc(1) * qc(4))];

        basis_c(:,3) = [2 * s * (qc(1) * qc(3) + qc(2) * qc(4)); ...
                        2 * s * (qc(2) * qc(3) - qc(1) * qc(4)); ...
                        1 - 2 * s * (qc(1)^2 + qc(2)^2)];

        % change of basis matrices
        C_Nc = C_cos_transf(eye(3), basis_c);
        
        % positions of the tyre connections
        pc1 = C_Nc * r1 + pcc;
        pc2 = C_Nc * r2 + pcc;
        pc3 = C_Nc * r3 + pcc;
        pc4 = C_Nc * r4 + pcc;
        
        y_return(i, 6) = pc1(2);
        y_return(i, 7) = pc2(2);
        y_return(i, 8) = pc3(2);
        y_return(i, 9) = pc4(2);
        
        % x and z coordinates of all legs
        y_return(i, 18) = pc1(1);
        y_return(i, 19) = pc2(1);
        y_return(i, 20) = pc3(1);
        y_return(i, 21) = pc4(1);

        y_return(i, 22) = pc1(3);
        y_return(i, 23) = pc2(3);
        y_return(i, 24) = pc3(3);
        y_return(i, 25) = pc4(3);


    end
end
