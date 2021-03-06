function y_return = calculateCornersMBD(y, r1, r2, r3, r4)
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
        
        pcc_inv = pcc;
        pcc_inv(3) = pcc(2);
        pcc_inv(2) = pcc(3);

        pc1 = C_Nc * r1 + pcc_inv;
        pc2 = C_Nc * r2 + pcc_inv;
        pc3 = C_Nc * r3 + pcc_inv;
        pc4 = C_Nc * r4 + pcc_inv;

        pw1 = y(i,38:40);
        pw2 = y(i,41:43);
        pw3 = y(i,44:46);
        pw4 = y(i,47:49);
        
        pt1 = y(i,50:52);
        pt2 = y(i,53:55);
        pt3 = y(i,56:58);
        pt4 = y(i,59:61);
        
        y_return(i,8) = pc1(1);        
        y_return(i,11) = pc2(1);
        y_return(i,14) = pc3(1);
        y_return(i,17) = pc4(1); 
        
        y_return(i,10) = pc1(2);        
        y_return(i,13) = pc2(2);
        y_return(i,16) = pc3(2);
        y_return(i,19) = pc4(2); 
        
        y_return(i,9) = pc1(3);        
        y_return(i,12) = pc2(3);
        y_return(i,15) = pc3(3);
        y_return(i,18) = pc4(3);  
        
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