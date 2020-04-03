clear; clc; close all;

syms dx1 dx2 k10 k20 Q11 Q12 Q21 Q22 q11 q12 q21 q22

Q = [Q11, Q12; Q21, Q22];
q = [q11, q12; q21, q22];

K = Q*[dx1^2; dx2^2] + q*[dx1; dx2] + [k10; k20];

f = [-K(1)*dx1 + K(2)*dx2; -K(2)*dx2];

J_dx = jacobian(f, [dx1;dx2]);
J_x_first = J_dx * [-1, 0; 1 -1];
%%
% clear; clc; close all;

syms dx1 dx2 k10 k20 Q11 Q12 Q21 Q22 q11 q12 q21 q22 x1 x2 l1 l2 m1 m2 g k1 k2

dx1 = l1 - x1;
dx2 = l2 -x2 +x1;

Q = [Q11, Q12; Q21, Q22];
q = [q11, q12; q21, q22];

K = Q*[dx1^2; dx2^2] + q*[dx1; dx2] + [k10; k20];

rhs = -m1*g+k1*l1-k2*l2;
f = [-K(1)+K(2), -K(2); -K(2), -K(2)] * [x1; x2] - rhs;
f1 = [-K(1)*dx1 + K(2)*dx2; -K(2)*dx2];

J_x = jacobian(f, [x1;x2]);
J_x1 = jacobian(f1, [x1;x2]);

simplify(J_x - J_x1)
simplify(J_x - J_x_first)
simplify(J_x_first - J_x1)