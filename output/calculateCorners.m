function y_return = calculateCorners(y, r1, r2, r3, r4)
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
%         pc1 = C_Nc * r1 + pcc;
%         pc2 = C_Nc * r2 + pcc;
%         pc3 = C_Nc * r3 + pcc;
%         pc4 = C_Nc * r4 + pcc;
        angles = zeros(3,1);
        yaw = y(i,33);
        pitch = y(i,32);
        roll = y(i,31);
        rotMat = smallAngleApproxMatrix(yaw, pitch, roll);
        
        pc1 = rotMat * r1 + pcc;
        pc2 = rotMat * r2 + pcc;
        pc3 = rotMat * r3 + pcc;
        pc4 = rotMat * r4 + pcc;

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