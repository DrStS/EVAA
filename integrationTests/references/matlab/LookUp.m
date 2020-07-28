% 
% script to create lookup tables
% - nonlinear stiffness: stiffness depends on length
%       displacement - force
% - nonlinear damping: stiffness depends on velocity
%       velocity - force

close all;
clear;
clc;

size_grid = 1000;

%params for stiffness lookup table:
stiff_name='U_Fk_lookup2.dat';
make_lookup_stiff=1;
ak = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3];
b=0;
c=0;
l_min = 0.05;
l_max = 0.8;
L = 0.3; %initial length of spring

%params for damping lookup table:
damp_name='V_Fd_lookup2.dat';
d=0.01;
make_lookup_damp=1;
ad = d* ak;
b_d=d*b;
c_d=d*c;
v_min=-1;
v_max=1;


if make_lookup_stiff
    dl = (l_max-l_min)/(size_grid-1);
    X = zeros(size_grid, 1);
    k_grid = zeros(size_grid,8);
    k = @(l,a)(c*l*l + b*l + a);
    for i = 0:size_grid-1
        X(i+1) = l_min+i*dl;
    end
    for i = 1:8 
        for j = 1:size_grid
            k_grid(j,i)= k(X(j), ak(i));
        end
    end
    u=X-L;
    fk=k_grid.*u;
    writematrix([u,fk],stiff_name,'Delimiter',',');
end


if make_lookup_damp
    dv = (v_max - v_min)/(size_grid-1);
    V = zeros(size_grid,1);
    d_grid = zeros(size_grid,8);
    d = @(l,a)(c_d*l*l + b_d *l + a);

    for i = 0:size_grid-1
        V(i+1) = v_min+i*dv;
    end
    for i = 1:8 
        for j = 1:size_grid
            d_grid(j,i)= d(V(j), ad(i));
        end
    end
    fd = d_grid.*V;
    writematrix([V,fd],damp_name,'Delimiter',',');   
end