function [q] = calculateQuaternions(angles, numIterations)
    q = zeros(numIterations, 4);
    for i=1:numIterations
        X = angles(i,1);        
        Y = angles(i,2);        
        Z = angles(i,3);    
        

        q(i,1) = sin(X/2)*cos(Y/2)*cos(Z/2) - cos(X/2)*sin(Y/2)*sin(Z/2);
        q(i,2) = cos(X/2)*sin(Y/2)*cos(Z/2) + sin(X/2)*cos(Y/2)*sin(Z/2);
        q(i,3) = cos(X/2)*cos(Y/2)*sin(Z/2) - sin(X/2)*sin(Y/2)*cos(Z/2);
        q(i,4) = cos(X/2)*cos(Y/2)*cos(Z/2) + sin(X/2)*sin(Y/2)*sin(Z/2);
    end