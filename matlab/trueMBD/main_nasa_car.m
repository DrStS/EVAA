function [t, y_return, y] = main_nasa_car(r1, r2, r3, r4, mass, mass_wheel, mass_tyre, Ic, qc, ... 
        lower_spring_length, upper_spring_length, initial_lower_spring_length, initial_upper_spring_length, ...
        lower_spring_stiffness, upper_spring_stiffness,  lower_spring_damping, upper_spring_damping, lower_rotational_stiffness, upper_rotational_stiffness, ...
        vc, vw1, vw2, vw3, vw4, vt1, vt2, vt3, vt4, wc, FC, FW1, FW2, FW3, FW4, FT1, FT2, FT3, FT4, FR1, FR2, FR3, FR4,...
        num_iter, delta_t, solver)
   
    %% INITIALLIZATION
    pcc = [0; 0; 0];                    % initial position of the center of the car at the origin
    
%     % positions of the wheels in global coordinates
%     pw1 = [r1(1);  r1(2)-initial_upper_spring_length(1); r1(3)  ];
%     pw2 = [r2(1);  r2(2)-initial_upper_spring_length(2); r2(3)  ];
%     pw3 = [r3(1);  r3(2)-initial_upper_spring_length(3); r3(3)  ];
%     pw4 = [r4(1);  r4(2)-initial_upper_spring_length(4); r4(3)  ];
%     
%     % position of the tyres in global coordinates
%     pt1 = [r1(1);  pw1(2)-initial_lower_spring_length(1)    ; r1(3)  ];
%     pt2 = [r2(1);  pw2(2)-initial_lower_spring_length(2)    ; r2(3)  ];
%     pt3 = [r3(1);  pw3(2)-initial_lower_spring_length(3)    ; r3(3)  ];
%     pt4 = [r4(1);  pw4(2)-initial_lower_spring_length(4)    ; r4(3)  ];
    [qc, pw1, pw2, pw3, pw4, pt1, pt2, pt3, pt4] = get_initial_positions(qc, r1, r2, r3, r4, pcc, initial_upper_spring_length, initial_lower_spring_length);


    % computation of ri_tilda, ro_tilda
    r1_tilda = get_tilda(r1);
    r2_tilda = get_tilda(r2);
    r3_tilda = get_tilda(r3);
    r4_tilda = get_tilda(r4);

    % initialise solution vector
    x_vector = [wc; ...     % 3 1:3 
                vc; ...     % 3 4:6
                vw1; ...    % 3 7:9
                vw2; ...    % 3 10:12
                vw3; ...    % 3 13:15
                vw4; ...    % 3 16:18
                vt1; ...    % 3 19:21
                vt2; ...    % 3 22:24
                vt3; ...    % 3 25:27
                vt4; ...    % 3 28:30
                qc; ...     % 4 31:34
                pcc; ...    % 3 35:37
                pw1; ...    % 3 38:40
                pw2; ...    % 3 41:43
                pw3; ...    % 3 44:46
                pw4; ...    % 3 47:49
                pt1; ...    % 3 50:52
                pt2; ...    % 3 53:55
                pt3; ...    % 3 56:58
                pt4];       % 3 59:61
    
    % create the reduced system matrix (probably best to store directly
    % its inverse, maybe split it into one vector of the inverse diagonal
    % elements and one matrix with the inverse of Ic)
    A = diag([1,1,1, mass, mass, mass, ...
        mass_wheel(1), mass_wheel(1), mass_wheel(1),...
        mass_wheel(2), mass_wheel(2), mass_wheel(2),...
        mass_wheel(3), mass_wheel(3), mass_wheel(3),...
        mass_wheel(4), mass_wheel(4), mass_wheel(4),...
        mass_tyre(1), mass_tyre(1), mass_tyre(1),...
        mass_tyre(2), mass_tyre(2), mass_tyre(2),...
        mass_tyre(3), mass_tyre(3), mass_tyre(3),...
        mass_tyre(4), mass_tyre(4), mass_tyre(4)]);  
    A(1:3, 1:3) = Ic;

    % create the structure with parameters used in functions' computation
    aux_vals = struct('r1', r1, ...
                      'r2', r2, ...
                      'r3', r3, ...
                      'r4', r4, ...
                      'r1_tilda', r1_tilda, ...
                      'r2_tilda', r2_tilda, ...
                      'r3_tilda', r3_tilda, ...
                      'r4_tilda', r4_tilda, ...
                      'FC', [0;FC;0], ...
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
                      'lower_spring_damping1', lower_spring_damping{1}, ...
                      'lower_spring_damping2', lower_spring_damping{2}, ...
                      'lower_spring_damping3', lower_spring_damping{3}, ...
                      'lower_spring_damping4', lower_spring_damping{4}, ...
                      'upper_spring_damping1', upper_spring_damping{1}, ...
                      'upper_spring_damping2', upper_spring_damping{2}, ...
                      'upper_spring_damping3', upper_spring_damping{3}, ...
                      'upper_spring_damping4', upper_spring_damping{4}, ...
                      'lower_rotational_stiffness', lower_rotational_stiffness, ...
                      'upper_rotational_stiffness', upper_rotational_stiffness, ...
                      'A', A);

    %% SOLVE PHASE              
    t = 0:delta_t:num_iter*delta_t;         

    f = @(t,x) compute_f3D_reduced(x,t,aux_vals);

    tic;
    [t,y] = solver(f, t, x_vector);
    solvetime=toc;
    timestr=['It took ', num2str(floor(solvetime)), 'sec ', num2str(round(mod(solvetime, 1)*1000)),'ms to solve the system!'];
    disp(timestr)
    disp(y(end,:))
    %% ONLY FOR VISUALISATION (NO CPP)
    y_return = zeros(size(y,1), 43);  % components: qc(1:4), pcc(5:7), 
                                      % pc1(8:10), pc2(11:13), pc3(14:16), pc4(17:19) 
                                      % pw1(20:22), pw2(23:25), pw3(26:28), pw4(29:31) 
                                      % pt1(32:34), pt2(35:37), pt3(38:40), pt4(41:43) 

    y_return(:,1:7) = y(:,31:37);
    % retrieve positions from the angles 
    for i = 1:size(y,1)
        qc = y_return(i,1:4)';
        pcc = y_return(i,5:7)'; 
       
        basis_c = get_basis(qc);

        % change of basis matrices
        C_Nc = C_cos_transf(eye(3), basis_c);
       
        % positions of the tyre connections
        pc1 = C_Nc * r1 + pcc;
        pc2 = C_Nc * r2 + pcc;
        pc3 = C_Nc * r3 + pcc;
        pc4 = C_Nc * r4 + pcc;

        pw1 = y(i,38:40);
        pw2 = y(i,41:43);
        pw3 = y(i,44:46);
        pw4 = y(i,47:49);
        
        pt1 = y(i,50:52);
        pt2 = y(i,53:55);
        pt3 = y(i,56:58);
        pt4 = y(i,59:61);
        
        y_return(i,8:10) = pc1;        
        y_return(i,11:13) = pc2;
        y_return(i,14:16) = pc3;
        y_return(i,17:19) = pc4;  
        
        y_return(i,20:22) = pw1;
        y_return(i,23:25) = pw2;
        y_return(i,26:28) = pw3;
        y_return(i,29:31) = pw4;
        
        y_return(i,32:34) = pt1;
        y_return(i,35:37) = pt2;
        y_return(i,38:40) = pt3;
        y_return(i,41:43) = pt4;

    end
end
