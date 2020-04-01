%clear; clc; close all;

%%
a = 1;
b = 2;
c = 3;
f1 = 1;
f2 = 0;
L1 = 1;
L2 = 1;

% parameters
tol = 1e-12;
x_init = [0; 0];
iter = 0;
err = [];

%% init declare variables

syms l1 u1 u2 l2 
k = @(a,b,c,l)(c*l*l + b*l + a);
l1 = u1 - u2 + L1;
l2 = u2 + L2;
k1 = k(a,b,c,l1);
k2 = k(a,b,c,l2);

%% Internal forces
syms p1 p2
p1 = k1*u1 - k2*u2;
p2 = -k1*u1+(k1+k2)*u2;

%% Jacobian
syms K k11 k12 k21 k22 K_symb
k11 = diff(p1, u1);
k12 = diff(p1, u2);
k21 = diff(p2, u1);
k22 = diff(p2, u2);

K_symb = [k11, k12; 
    k21, k22];

r_symb = [p1 - f1; 
     p2 - f2];
 
%% Newton iteration
u1 = 0;
u2 = 0;
u = [u1; u2];
Delta = [u1; u2];

while 1
    iter = iter + 1;
    
   K = eval(K_symb);
   r = eval(r_symb);
   Delta = -K\r;
   u = Delta + u;
   u1 = u(1);
   u2 = u(2);
   
   err = [err, norm(r)];
   if (err(length(err)) < tol)
        break;
   end
   
end
iter
%order = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))))
disp(u)
disp(err(length(err)))