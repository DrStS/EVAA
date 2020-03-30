clear; clc; close all;

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

%% interpolation
global dl k_grid l_min
% evaluation density
dl = 0.1;
% calc min and max length to evaluate (just for now)
l_min = min(0.3 * L1,0.3 * L2);
l_max = max(1.7 * L1,1.7 * L2);
% grid size
k_d = round((l_max - l_min)/dl);
% grid value allocation
k_grid = zeros(k_d+1,1);
% k response function
k = @(l)(c*l*l + b*l + a);
% fill in grid values
for i = 0:k_d
        k_grid(i+1)= k(l_min+i*dl);
end

%% init declare variables

syms l1 u1 u2 l2 ;
l1 = u1 - u2 + L1;
l2 = u2 + L2;

%% Internal forces
syms p1 p2 k1 k2;
p1 = k1*u1 - k2*u2;
p2 = -k1*u1+(k1+k2)*u2;

%% Jacobian
syms K k11 k12 k21 k22 K_symb
% actuall derivative but sym cant be used as we get part of derivative out
% my function
%k11 = diff(p1, u1)+u1*deriv_dk(l1);
%k12 = diff(p1, u2)-u1*deriv_dk(l1)-u2*deriv_dk(l2);
%k21 = diff(p2, u1)+(u2-u1)*deriv_dk(l1);
%k22 = diff(p2, u2)-(u2-u1)*deriv_dk(l1)+u2*deriv_dk(l2);
k11 = diff(p1, u1);
k12 = diff(p1, u2);
k21 = diff(p2, u1);
k22 = diff(p2, u2);

K_symb_1 = [k11, k12; 
    k21, k22];

r_symb = [p1 - f1; 
     p2 - f2];

 
%% Newton iteration
u1 = 0;
u2 = 0;
u = [u1; u2];
Delta = [u1; u2];

% inital values
l_1 = eval(l1);
l_2 = eval(l2);
k1 = interpolate(l_1);
k2 = interpolate(l_2);

while 1
   iter = iter + 1;
   K = eval(K_symb_1)+[u(1),-u(1);u(2)-u(1),u(1)-u(2)]*deriv_dk(l_1)+ [0 -u(2); 0 u(2)]*deriv_dk(l_2);
   r = eval(r_symb);
   Delta = -K\r;
   u = Delta + u;
   % update values
   u1 = u(1);
   u2 = u(2);
   l_1 = eval(l1);
   l_2 = eval(l2);
   k1 = interpolate(l_1);
   k2 = interpolate(l_2);
   
   err = [err, norm(r)];
   if (err(length(err)) < tol)
        break;
   end
end
iter
order = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))))
disp(u)
disp(err(length(err)))
ref = [ 0.2788; 0.1394]

% interpolates k_grid to a certain value l
function k = interpolate(l)
    global k_grid l_min dl;
    % matlab need a +1 
    % we dont need it in c++ tho
    idx = (l-l_min)/dl + 1;
    if (floor(idx) == idx)
        k = k_grid(idx);
    else
        i_low = floor(idx);
        low = k_grid(i_low);
        top = k_grid(ceil(idx));
        pr = (l-(l_min+dl*(i_low-1)))/dl; % percentage from left to right
        k = (1-pr)*low+pr*top;
    end
end

% calculates the derivative at a certain position of the grid
function dk = deriv_dk(l)
    global k_grid l_min dl;
    idx = (l-l_min)/dl + 1;
    if (floor(idx) == idx)
        low = k_grid(idx-1);
        top = k_grid(idx+1);
        dk = (top-low)/(2*dl);
    else
        low = k_grid(floor(idx));
        top = k_grid(ceil(idx));
        dk = (top-low)/dl;
    end
end