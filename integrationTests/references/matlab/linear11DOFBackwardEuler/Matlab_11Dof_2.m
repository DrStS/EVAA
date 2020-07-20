%%
% Copyright &copy; 2020, Nicola Zimmermann, Munich \n
% All rights reserved. \n
% 
% This file is part of EVAA.
% 
% EVAA is free software: you can redistribute it and/or modify \n
% it under the terms of the GNU General Public  License as published by \n
% the Free Software Foundation, either version 3 of the License, or \n
% (at your option) any later version. \n
% 
% EVAA is distributed in the hope that it will be useful, \n
% but WITHOUT ANY WARRANTY; without even the implied warranty of \n
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the \n
% GNU General Public License for more details. \n
% 
% You should have received a copy of the GNU General Public License \n
% along with EVAA.  If not, see <a
% href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses/</a>.
% 
% Additional permission under GNU GPL version 3 section 7
% 
% If you modify this Program, or any covered work, by linking or combining it
% with Intel Math Kernel Libraries(MKL) (or a modified version of that
% library), containing parts covered by the terms of the license of the MKL,
% the licensors of this Program grant you additional permission to convey the
% resulting work.
% 
% DESCRIPTION
% This file is used to compute displacements/ rotations
% of the 11-Dof 2-track-model (Eulerian Dofs)
% with a backward Euler time integration scheme 
% for nonlinear spring stiffness (depends on length)
% stiffness values (and derivatives) are obtained from lookup table (precomputed)
% with linear or spline interpolation
% stiffness set to constant

%%
clear; 
% % clc; 
close all;
%% make lookup table
ak = [19.32e3; 260e3; 19.32e3; 260e3; 13.12e3; 260e3; 13.12e3; 260e3];
% b = 1.59;
% c = 1.2;
b=0;
c=0;

global k_grid X size_grid;
l_min = 0.05;
l_max = 0.8;
size_grid = 10;
X = zeros(size_grid, 1);
k_grid = zeros(8*size_grid,1);
dl = (l_max-l_min)/(size_grid-1);
k = @(l,a)(c*l*l + b*l + a);
for i = 0:size_grid-1
    X(i+1) = l_min+i*dl;
end
for i = 1:8 
    for j = 1:size_grid
        k_grid((i-1)*size_grid + j)= k(X(j), ak(i));
    end
end
 
global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 

k_spline1 = spline(X,k_grid(1:size_grid));
k_der1 = fnder(k_spline1,1);
k_spline2 = spline(X,k_grid(size_grid+1:2*size_grid));
k_der2 = fnder(k_spline2,1);
k_spline3 = spline(X,k_grid(2*size_grid+1:3*size_grid));
k_der3 = fnder(k_spline3,1);
k_spline4 = spline(X,k_grid(3*size_grid+1:4*size_grid));
k_der4 = fnder(k_spline4,1);
k_spline5 = spline(X,k_grid(4*size_grid+1:5*size_grid));
k_der5 = fnder(k_spline5,1);
k_spline6 = spline(X,k_grid(5*size_grid+1:6*size_grid));
k_der6 = fnder(k_spline6,1);
k_spline7 = spline(X,k_grid(6*size_grid+1:7*size_grid));
k_der7 = fnder(k_spline7,1);
k_spline8 = spline(X,k_grid(7*size_grid+1:8*size_grid));
k_der8 = fnder(k_spline8,1);


%% parameters
global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr
global l1 l2 l3 l4 l5 l6 l7 l8
global p f
global K dKdxx



num_timesteps = 1000;
delta_t = 1/num_timesteps; 

tol = 1e-7;

l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=1.6916;
l_lat_fr=1.6916;
l_lat_rl=1.68;
l_lat_rr=1.68;
mass_Body=1936;
I_body_xx=640;
I_body_yy=4800;
mass_wheel_fl=72.5;
mass_tyre_fl=0;
mass_wheel_fr=72.5;
mass_tyre_fr=0;
mass_wheel_rl=67.5;
mass_tyre_rl=0;
mass_wheel_rr=67.5;
mass_tyre_rr=0;

u_n = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

f =[1.1e3; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

L = 0.3;
%interpolation type in lookup table (linear/spline)
type='spline';

%% symbolic derivative

% syms l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr
syms k1 k2 k3 k4 k5 k6 k7 k8
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
syms d1 d2 d3 d4 d5 d6 d7 d8

K = [k1+k3+k5+k7, -k1*l_lat_fl+k3*l_lat_fr-k5*l_lat_rl+k7*l_lat_rr, -k1*l_long_fl-k3*l_long_fr+k5*l_long_rl+k7*l_long_rr,  -k1, 0, -k3, 0, -k5, 0, -k7, 0;
   0, l_lat_fl*l_lat_fl*k1+l_lat_fr*l_lat_fr*k3+l_lat_rl*l_lat_rl*k5+l_lat_rr*l_lat_rr*k7, +l_long_fl*l_lat_fl*k1-l_lat_fr*l_long_fr*k3-l_long_rl*l_lat_rl*k5+l_long_rr*l_lat_rr*k7, +l_lat_fl*k1, 0, -l_lat_fr*k3, 0, +l_lat_rl*k5, 0, -l_lat_rr*k7, 0;
   0, 0, l_long_fl*l_long_fl*k1+l_long_fr*l_long_fr*k3+l_long_rl*l_long_rl*k5+l_long_rr*l_long_rr*k7, l_long_fl*k1, 0, l_long_fr*k3, 0, -l_long_rl*k5, 0, -l_long_rr*k7, 0; 
   0, 0, 0, k1+k2, -k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k3+k4, -k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k5+k6, -k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k7+k8, -k8;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k8];
K=K+K.'-diag(diag(K));
   
x=[x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11];
p=K*x;

l1=diff(p(1),k1)+L;
l2=diff(p(4),k2)+L;
l3=diff(p(1),k3)+L;
l4=diff(p(6),k4)+L;
l5=diff(p(1),k5)+L;
l6=diff(p(8),k6)+L;
l7=diff(p(1),k7)+L;
l8=diff(p(10),k8)+L;

k1=d1*l1;
k2=d2*l2;
k3=d3*l3;
k4=d4*l4;
k5=d5*l5;
k6=d6*l6;
k7=d7*l7;
k8=d8*l8;

dKdxx=diff(eval(K),x1)*x1+diff(eval(K),x2)*x2...
    + diff(eval(K),x3)*x3+diff(eval(K),x4)*x4...
    + diff(eval(K),x5)*x5+diff(eval(K),x6)*x6...
    + diff(eval(K),x7)*x7+diff(eval(K),x8)*x8...
    + diff(eval(K),x9)*x9+diff(eval(K),x10)*x10...
    + diff(eval(K),x11)*x11;


%%

M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
M_div_h2 = M / (delta_t * delta_t);

residfunc = @(y_curr,y1,y2,K,rhs)( ( M_div_h2 + K ) * y_curr - 2 * M_div_h2 * y1 + M_div_h2 * y2 - rhs);


%% Time Loop

t = 0:delta_t:(num_timesteps)*delta_t;

u_n_m_1 = u_n;
u_n_p_1=u_n;

sol_mat = zeros(length(t),11);
newt_order_vec = zeros(length(t),1);
err_vec = zeros(length(t),1);
condition = zeros(length(t),1);
iterat_mat=zeros(length(t),1);

%init values
x1 = u_n(1);
x2 = u_n(2);
x3 = u_n(3);
x4 = u_n(4);
x5 = u_n(5);
x6 = u_n(6);
x7 = u_n(7);
x8 = u_n(8);
x9 = u_n(9);
x10 = u_n(10);
x11 = u_n(11);

sol_mat(1,:)=u_n;

for i = 2: length(t)
    K1=kfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type);
    u_n_p_1 = ( M_div_h2 + K1 )\ (2 * M_div_h2 * u_n - M_div_h2 * u_n_m_1 + f);
    x1 = u_n_p_1(1);
    x2 = u_n_p_1(2);
    x3 = u_n_p_1(3);
    x4 = u_n_p_1(4);
    x5 = u_n_p_1(5);
    x6 = u_n_p_1(6);
    x7 = u_n_p_1(7);
    x8 = u_n_p_1(8);
    x9 = u_n_p_1(9);
    x10 = u_n_p_1(10);
    x11 = u_n_p_1(11);
    

    K1=kfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type);
    r=residfunc(u_n_p_1,u_n,u_n_m_1,K1,f);
    
    err = [];
    delta = [];
        
    for iter=1:10

        K2 = kderivfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type);        
        J = M_div_h2 + K1 + K2;
%         r=residfunc(u_n,u_n_p_1,u_n_m_1,K1,f);       
        Delta = -J\r;
        u_n_p_1 = Delta + u_n_p_1; 

        x1 = u_n_p_1(1);
        x2 = u_n_p_1(2);
        x3 = u_n_p_1(3);
        x4 = u_n_p_1(4);
        x5 = u_n_p_1(5);
        x6 = u_n_p_1(6);
        x7 = u_n_p_1(7);
        x8 = u_n_p_1(8);
        x9 = u_n_p_1(9);
        x10 = u_n_p_1(10);
        x11 = u_n_p_1(11);

        K1=kfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type); 
        r=residfunc(u_n_p_1,u_n,u_n_m_1,K1,f);
        
        err(iter)=norm(r);
        if (err(iter) < tol || norm(J\r) > norm(Delta))
            err_vec(i)=norm(r);
            break;
        end       
    end
    
    if iter>2
        newt_order_vec(i) = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))));
    end
    iterat_mat(i) = iter;
    condition(i) =  cond(J);
    
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
    sol_mat(i,:) = u_n_p_1;
end

disp(sol_mat(1001,1:3))

figure;
subplot(1,3,1);
plot(t,sol_mat(:,1)); grid on; legend('z_{CG}'); 
subplot(1,3,2); 
plot(t,sol_mat(:,2)); grid on; legend('r_x');
subplot(1,3,3);
plot(t,sol_mat(:,3)); grid on; legend('r_y');

figure; plot(err_vec); title('error'); grid on;
figure; plot(iterat_mat); title('newton iterations'); grid on;

%% write to textfile
writematrix([t',sol_mat],'Matlab_11Dof_2.dat','Delimiter',',');

%%
function K1 = kfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type)      
    
    global X k_grid size_grid
    global K
    global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr 
    global l1 l2 l3 l4 l5 l6 l7 l8
    global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
    
    switch type
        case 'linear'
            k1=interp1(X,k_grid(1:size_grid),eval(l1),'linear',k_grid(end));
            k2=interp1(X,k_grid(size_grid+1:2*size_grid),eval(l2),'linear',k_grid(end));          
            k3=interp1(X,k_grid(2*size_grid+1:3*size_grid),eval(l3),'linear',k_grid(end));
            k4=interp1(X,k_grid(3*size_grid+1:4*size_grid),eval(l4),'linear',k_grid(end));
            k5=interp1(X,k_grid(4*size_grid+1:5*size_grid),eval(l5),'linear',k_grid(end));
            k6=interp1(X,k_grid(5*size_grid+1:6*size_grid),eval(l6),'linear',k_grid(end));          
            k7=interp1(X,k_grid(6*size_grid+1:7*size_grid),eval(l7),'linear',k_grid(end));
            k8=interp1(X,k_grid(7*size_grid+1:end),eval(l8),'linear',k_grid(end)); 
            
        case 'spline'
            k1=ppval(k_spline1,eval(l1));
            k2=ppval(k_spline2,eval(l2));
            k3=ppval(k_spline3,eval(l3));
            k4=ppval(k_spline4,eval(l4));
            k5=ppval(k_spline5,eval(l5));
            k6=ppval(k_spline6,eval(l6));
            k7=ppval(k_spline7,eval(l7));
            k8=ppval(k_spline8,eval(l8));
    
    end   
    K1=eval(K);
end

function K2 = kderivfunc(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,type)      
    
    global X k_grid size_grid
    global dKdxx
    global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr 
    global l1 l2 l3 l4 l5 l6 l7 l8
    global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 
    
    switch type
        case 'linear'
            k1=interp1(X,k_grid(1:size_grid),eval(l1),'linear',k_grid(end));
            k2=interp1(X,k_grid(size_grid+1:2*size_grid),eval(l2),'linear',k_grid(end));          
            k3=interp1(X,k_grid(2*size_grid+1:3*size_grid),eval(l3),'linear',k_grid(end));
            k4=interp1(X,k_grid(3*size_grid+1:4*size_grid),eval(l4),'linear',k_grid(end));
            k5=interp1(X,k_grid(4*size_grid+1:5*size_grid),eval(l5),'linear',k_grid(end));
            k6=interp1(X,k_grid(5*size_grid+1:6*size_grid),eval(l6),'linear',k_grid(end));          
            k7=interp1(X,k_grid(6*size_grid+1:7*size_grid),eval(l7),'linear',k_grid(end));
            k8=interp1(X,k_grid(7*size_grid+1:end),eval(l8),'linear',k_grid(end)); 
            d1=kderivfun(k1,eval(l1),1);
            d2=kderivfun(k2,eval(l2),2);
            d3=kderivfun(k3,eval(l3),3);
            d4=kderivfun(k4,eval(l4),4);
            d5=kderivfun(k5,eval(l5),5);
            d6=kderivfun(k6,eval(l6),6);
            d7=kderivfun(k7,eval(l7),7);
            d8=kderivfun(k8,eval(l8),8);
            
        case 'spline'
            d1=ppval(k_der1,eval(l1));
            d2=ppval(k_der2,eval(l2));
            d3=ppval(k_der3,eval(l3));
            d4=ppval(k_der4,eval(l4));
            d5=ppval(k_der5,eval(l5));
            d6=ppval(k_der6,eval(l6));
            d7=ppval(k_der7,eval(l7));
            d8=ppval(k_der8,eval(l8));
    
    end   
    K2=eval(dKdxx); 
end

function dk_du = kderivfun(k,l,idx)
    global X k_grid size_grid 

    k_grid_tmp=k_grid((idx-1)*size_grid+1:idx*size_grid);
    knext=interp1(X,k_grid_tmp,l,'next',k_grid(1));

    if k==knext && k_grid_tmp(1)==k_grid_tmp(2)
        dk_du=0;
    elseif k==knext && k_grid_tmp(1)~= k_grid_tmp(2) 
        knext=k_grid_tmp(find(k_grid_tmp==k)+1);
        kprev=k_grid_tmp(find(k_grid_tmp==k)-1);
        lnext=X(find(k_grid_tmp==knext));
        lprev=X(find(k_grid_tmp==kprev));
        dk_du=(knext-kprev)/(lnext-lprev); 
    else
        lnext=X(k_grid_tmp==knext);
        dk_du=(knext-k)/(lnext-l);
    end
end

