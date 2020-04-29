function [basis] = get_basis(q)
% quaternion q yields basis (calculates basically the matrix rotation
% the euclidian basis to the local basis)

        %get local basis vectors using unit quaternion rotation 
        s = 1 / norm(q)^2;      %normalizer, only to remove numerical stuffy stuff

        basis = zeros(3);
        basis(:,1) = [1 - 2 * s *(q(2)^2 + q(3)^2); ...
                        2 * s * (q(1) * q(2) + q(3) * q(4)); ...
                        2 * s * (q(1) * q(3) - q(2) * q(4))];

        basis(:,2) = [2 * s * (q(1) * q(2) - q(3) * q(4)); ...
                        1 - 2 * s * (q(1)^2 + q(3)^2); ...
                        2 * s * (q(2) * q(3) + q(1) * q(4))];

        basis(:,3) = [2 * s * (q(1) * q(3) + q(2) * q(4)); ...
                        2 * s * (q(2) * q(3) - q(1) * q(4)); ...
                        1 - 2 * s * (q(1)^2 + q(2)^2)];
end

