function C = C_cos_transf(X, Y)
% This function defines the cosine transformation matrix, used in changing
% coordinates from a frame to another.
% The X and Y are matrices the represent basis in the X and, respectively,
% Y frame.
%   x = C * y

    C = zeros(3);
    for i = 1 : 3
        for j = 1 : 3
            C(i,j) = dot(X(:,i),Y(:,j))/(norm(X(:,i))*norm(Y(:,j)));
        end
    end
end