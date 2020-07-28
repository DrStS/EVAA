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
% for nonlinear spring and damping
% stiffness & damping values (and derivatives) are obtained from lookup table
% lookup table gives displ/vel & force pairs (precomputed)

%%
clear; 
% % clc; 
% close all;

%% parameters
% global l1 l2 l3 l4 l5 l6 l7 l8
% global K dKdxx
% global D dDdxv

num_timesteps = 1000;
delta_t = 1/num_timesteps; 

tol = 1e-7;

k_stab_f=20e3;
k_stab_r=4e3;
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

syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
syms v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11
syms k1 k2 k3 k4 k5 k6 k7 k8
syms d_k1 d_k2 d_k3 d_k4 d_k5 d_k6 d_k7 d_k8
syms d1 d2 d3 d4 d5 d6 d7 d8
syms d_d1 d_d2 d_d3 d_d4 d_d5 d_d6 d_d7 d_d8


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

K_stab=[0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 -k_stab_f 0 k_stab_f 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 -k_stab_f 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 -k_stab_r 0 k_stab_r 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 -k_stab_r 0;
    0 0 0 0 0 0 0 0 0 0 0];

K=K+K_stab;

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

k1=d_k1*l1;
k2=d_k2*l2;
k3=d_k3*l3;
k4=d_k4*l4;
k5=d_k5*l5;
k6=d_k6*l6;
k7=d_k7*l7;
k8=d_k8*l8;

dKdxx=diff(eval(K),x1)*x1+diff(eval(K),x2)*x2...
    + diff(eval(K),x3)*x3+diff(eval(K),x4)*x4...
    + diff(eval(K),x5)*x5+diff(eval(K),x6)*x6...
    + diff(eval(K),x7)*x7+diff(eval(K),x8)*x8...
    + diff(eval(K),x9)*x9+diff(eval(K),x10)*x10...
    + diff(eval(K),x11)*x11;

%damping
D = [d1+d3+d5+d7, -d1*l_lat_fl+d3*l_lat_fr-d5*l_lat_rl+d7*l_lat_rr, -d1*l_long_fl-d3*l_long_fr+d5*l_long_rl+d7*l_long_rr,  -d1, 0, -d3, 0, -d5, 0, -d7, 0;
   0, l_lat_fl*l_lat_fl*d1+l_lat_fr*l_lat_fr*d3+l_lat_rl*l_lat_rl*d5+l_lat_rr*l_lat_rr*d7, +l_long_fl*l_lat_fl*d1-l_lat_fr*l_long_fr*d3-l_long_rl*l_lat_rl*d5+l_long_rr*l_lat_rr*d7, +l_lat_fl*d1, 0, -l_lat_fr*d3, 0, +l_lat_rl*d5, 0, -l_lat_rr*d7, 0;
   0, 0, l_long_fl*l_long_fl*d1+l_long_fr*l_long_fr*d3+l_long_rl*l_long_rl*d5+l_long_rr*l_long_rr*d7, l_long_fl*d1, 0, l_long_fr*d3, 0, -l_long_rl*d5, 0, -l_long_rr*d7, 0; 
   0, 0, 0, d1+d2, -d2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, d2, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, d3+d4, -d4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, d4, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, d5+d6, -d6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, d6, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, d7+d8, -d8;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, d8];
D=D+D.'-diag(diag(D));

d1=d_d1*l1;
d2=d_d2*l2;
d3=d_d3*l3;
d4=d_d4*l4;
d5=d_d5*l5;
d6=d_d6*l6;
d7=d_d7*l7;
d8=d_d8*l8;

dDdxv=diff(eval(D),x1)*v1+diff(eval(D),x2)*v2...
    + diff(eval(D),x3)*v3+diff(eval(D),x4)*v4...
    + diff(eval(D),x5)*v5+diff(eval(D),x6)*v6...
    + diff(eval(D),x7)*v7+diff(eval(D),x8)*v8...
    + diff(eval(D),x9)*v9+diff(eval(D),x10)*v10...
    + diff(eval(D),x11)*v11;

K_fct=matlabFunction(K,'File','fct_K.m');
dK_fct=matlabFunction(dKdxx,'File','fct_dK.m');
D_fct=matlabFunction(D,'File','fct_D.m');
dD_fct=matlabFunction(dDdxv,'File','fct_dD.m');

clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11
clear v1 v2 v3 v4 v5 v6 v7 v8 v9 v10 v11
clear k1 k2 k3 k4 k5 k6 k7 k8
clear d_k1 d_k2 d_k3 d_k4 d_k5 d_k6 d_k7 d_k8
clear d1 d2 d3 d4 d5 d6 d7 d8
clear d_d1 d_d2 d_d3 d_d4 d_d5 d_d6 d_d7 d_d8
clear K dKdxx D dDdxv

%%

M = diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
M_div_h2 = M / (delta_t * delta_t);

residfunc = @(y_p1,y,y_m1,K,D,f)( ( M_div_h2 + D/delta_t + K ) * y_p1 - (2 * M_div_h2 + D/delta_t)* y + M_div_h2 * y_m1 - f);


%% Time Loop

t = 0:delta_t:(num_timesteps)*delta_t;

u_n_m_1 = u_n;
u_n_p_1=u_n;

sol_mat = zeros(length(t),11);
newt_order_vec = zeros(length(t),1);
err_vec = zeros(length(t),1);
condition = zeros(length(t),1);
iterat_mat=zeros(length(t),1);

sol_mat(1,:)=u_n;

for i = 2: length(t)
    K1=kfunc(u_n,type);
    D1=dfunc(u_n_p_1,u_n,delta_t,type);
    u_n_p_1 = ( M_div_h2 + D1/delta_t + K1 )\ ((2 * M_div_h2 + D1/delta_t) * u_n - M_div_h2 * u_n_m_1 + f);
    

    K1=kfunc(u_n_p_1,type);
    D1=dfunc(u_n_p_1,u_n,delta_t,type);
    r=residfunc(u_n_p_1,u_n,u_n_m_1,K1,D1,f);
    
    err = [];
    delta = [];
        
    for iter=1:10

        K2 = kderivfunc(u_n_p_1,type);        
        D2 = dderivfunc(u_n_p_1,u_n,delta_t,type);
        J = M_div_h2 + K1 + K2 + D1/delta_t + D2;
%         r=residfunc(u_n,u_n_p_1,u_n_m_1,K1,f);       
        Delta = -J\r;
        u_n_p_1 = Delta + u_n_p_1; 

        K1=kfunc(u_n_p_1,type);
        D1=dfunc(u_n_p_1,u_n,delta_t,type);
        r=residfunc(u_n_p_1,u_n,u_n_m_1,K1,D1,f);
        
        err(iter)=norm(r);
        if (err(iter) < tol || norm(J\r) > norm(Delta))
            err_vec(i)=norm(r);
            break;
        end       
    end
    
    if iter>3
        newt_order_vec(i) = log(abs((err(end)-err(end-1))/(err(end-1)-err(end-2))))/log(abs((err(end-1)-err(end-2))/(err(end-2)-err(end-3))));
    end
    iterat_mat(i) = iter;
    condition(i) =  cond(J);
    
    u_n_m_1 = u_n;
    u_n = u_n_p_1;
    sol_mat(i,:) = u_n_p_1;
    
%     if sum([eval(l1),eval(l2),eval(l3),eval(l4),eval(l5),eval(l6),eval(l7),eval(l8)]>l_max)
%         disp('Warning: l larger than l_max in iteration'); disp(i);
%     end
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
writematrix([t',sol_mat],'Matlab_11Dof_3dn.dat','Delimiter',',');

%%
function K1 = kfunc(u,type)      

    l=dufunc(u);
    
    tmp=readmatrix('U_fk_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    u_vec=tmp(:,1);
    k_mat=tmp(:,2:end)./u_vec;
    l_vec=u_vec;
    clear tmp;
       
    k_spline1 = spline(l_vec,k_mat(:,1));
    k_spline2 = spline(l_vec,k_mat(:,2));
    k_spline3 = spline(l_vec,k_mat(:,3));
    k_spline4 = spline(l_vec,k_mat(:,4));
    k_spline5 = spline(l_vec,k_mat(:,5));
    k_spline6 = spline(l_vec,k_mat(:,6));
    k_spline7 = spline(l_vec,k_mat(:,7));
    k_spline8 = spline(l_vec,k_mat(:,8));

%     global X k_grid size_grid
%     global K
%     global l1 l2 l3 l4 l5 l6 l7 l8
%     global k_spline1 k_spline2 k_spline3 k_spline4 k_spline5 k_spline6 k_spline7 k_spline8;
    
    switch type
        case 'linear'
            k1=interp1(u_vec,k_mat(:,1),l(1),'linear',k_mat(end,1));
            k2=interp1(u_vec,k_mat(:,2),l(2),'linear',k_mat(end,2));
            k3=interp1(u_vec,k_mat(:,3),l(3),'linear',k_mat(end,3));
            k4=interp1(u_vec,k_mat(:,4),l(4),'linear',k_mat(end,4));
            k5=interp1(u_vec,k_mat(:,5),l(5),'linear',k_mat(end,5));
            k6=interp1(u_vec,k_mat(:,6),l(6),'linear',k_mat(end,6));
            k7=interp1(u_vec,k_mat(:,7),l(7),'linear',k_mat(end,7));
            k8=interp1(u_vec,k_mat(:,8),l(8),'linear',k_mat(end,8));
               
        case 'spline'
            k1=ppval(k_spline1,l(1));
            k2=ppval(k_spline2,l(2));
            k3=ppval(k_spline3,l(3));
            k4=ppval(k_spline4,l(4));
            k5=ppval(k_spline5,l(5));
            k6=ppval(k_spline6,l(6));
            k7=ppval(k_spline7,l(7));
            k8=ppval(k_spline8,l(8));
    
    end   
%     K1=eval(K);
    K1=fct_K(k1,k2,k3,k4,k5,k6,k7,k8);
end

function K2 = kderivfunc(u,type)      
    
    x1=u(1);
    x2=u(2);
    x3=u(3);
    x4=u(4);
    x5=u(5);
    x6=u(6);
    x7=u(7);
    x8=u(8);
    x9=u(9);
    x10=u(10);
    x11=u(11);
    
    l=dufunc(u);

    tmp=readmatrix('U_fk_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    u_vec=tmp(:,1);
    k_mat=tmp(:,2:end)./u_vec;
    l_vec=u_vec;
    clear tmp;
       
    k_spline1 = spline(l_vec,k_mat(:,1));
    k_spline2 = spline(l_vec,k_mat(:,2));
    k_spline3 = spline(l_vec,k_mat(:,3));
    k_spline4 = spline(l_vec,k_mat(:,4));
    k_spline5 = spline(l_vec,k_mat(:,5));
    k_spline6 = spline(l_vec,k_mat(:,6));
    k_spline7 = spline(l_vec,k_mat(:,7));
    k_spline8 = spline(l_vec,k_mat(:,8));
    k_der1=fnder(k_spline1,1);
    k_der2=fnder(k_spline2,1);
    k_der3=fnder(k_spline3,1);
    k_der4=fnder(k_spline4,1);
    k_der5=fnder(k_spline5,1);
    k_der6=fnder(k_spline6,1);
    k_der7=fnder(k_spline7,1);
    k_der8=fnder(k_spline8,1);

%     global X k_grid size_grid
%     global dKdxx
%     global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr 
%     global l1 l2 l3 l4 l5 l6 l7 l8
%     global k_der1 k_der2 k_der3 k_der4 k_der5 k_der6 k_der7 k_der8; 
    
    switch type
        case 'linear'
            k1=interp1(u_vec,k_mat(:,1),l(1),'linear',k_mat(end,1));
            k2=interp1(u_vec,k_mat(:,2),l(2),'linear',k_mat(end,2));
            k3=interp1(u_vec,k_mat(:,3),l(3),'linear',k_mat(end,3));
            k4=interp1(u_vec,k_mat(:,4),l(4),'linear',k_mat(end,4));
            k5=interp1(u_vec,k_mat(:,5),l(5),'linear',k_mat(end,5));
            k6=interp1(u_vec,k_mat(:,6),l(6),'linear',k_mat(end,6));
            k7=interp1(u_vec,k_mat(:,7),l(7),'linear',k_mat(end,7));
            k8=interp1(u_vec,k_mat(:,8),l(8),'linear',k_mat(end,8));
            d_k1=kderivfun(k1,l(1),1);
            d_k2=kderivfun(k2,l(2),2);
            d_k3=kderivfun(k3,l(3),3);
            d_k4=kderivfun(k4,l(4),4);
            d_k5=kderivfun(k5,l(5),5);
            d_k6=kderivfun(k6,l(6),6);
            d_k7=kderivfun(k7,l(7),7);
            d_k8=kderivfun(k8,l(8),8);
            
        case 'spline'
            d_k1=ppval(k_der1,l(1));
            d_k2=ppval(k_der2,l(2));
            d_k3=ppval(k_der3,l(3));
            d_k4=ppval(k_der4,l(4));
            d_k5=ppval(k_der5,l(5));
            d_k6=ppval(k_der6,l(6));
            d_k7=ppval(k_der7,l(7));
            d_k8=ppval(k_der8,l(8));
    
    end   
    K2=fct_dK(d_k1,d_k2,d_k3,d_k4,d_k5,d_k6,d_k7,d_k8,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11); 
end

function dk_du = kderivfun(k,l,idx)
    
    tmp=readmatrix('U_fk_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    u_vec=tmp(:,1);
    k_mat=tmp(:,2:end)./u_vec;
    clear tmp;

% global X k_grid size_grid 

    k_mat_tmp=k_mat(:,idx);
    knext=interp1(u_vec,k_mat_tmp,l,'next',k_mat_tmp(end));

    if k==knext && k_mat_tmp(1)==k_mat_tmp(2)
        dk_du=0;
    elseif k==knext && k_mat_tmp(1)~= k_mat_tmp(2) 
        knext=k_mat_tmp(find(k_mat_tmp==k)+1);
        kprev=k_mat_tmp(find(k_mat_tmp==k)-1);
        lnext=u_vec(find(k_mat_tmp==knext));
        lprev=u_vec(find(k_mat_tmp==kprev));
        dk_du=(knext-kprev)/(lnext-lprev); 
    else
        lnext=u_vec(k_mat_tmp==knext);
        dk_du=(knext-k)/(lnext-l);
    end
   
end

function D1 = dfunc(u_n_p_1,u_n,delta_t,type)      
    
    v=1/delta_t.*(u_n_p_1-u_n);
    
    dv=dufunc(v);
    
    tmp=readmatrix('V_Fd_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    v_vec=tmp(:,1);
    d_mat=tmp(:,2:end)./v_vec;
    clear tmp;
       
    d_spline1 = spline(v_vec,d_mat(:,1));
    d_spline2 = spline(v_vec,d_mat(:,2));
    d_spline3 = spline(v_vec,d_mat(:,3));
    d_spline4 = spline(v_vec,d_mat(:,4));
    d_spline5 = spline(v_vec,d_mat(:,5));
    d_spline6 = spline(v_vec,d_mat(:,6));
    d_spline7 = spline(v_vec,d_mat(:,7));
    d_spline8 = spline(v_vec,d_mat(:,8));

%     global X d_grid size_grid
%     global D
%     global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr 
%     global l1 l2 l3 l4 l5 l6 l7 l8
%     global d_spline1 d_spline2 d_spline3 d_spline4 d_spline5 d_spline6 d_spline7 d_spline8;
    
    switch type
        case 'linear'
            d1=interp1(v_vec,d_mat(:,1),dv(1),'linear',d_mat(end,1));
            d2=interp1(v_vec,d_mat(:,2),dv(2),'linear',d_mat(end,2));
            d3=interp1(v_vec,d_mat(:,3),dv(3),'linear',d_mat(end,3));
            d4=interp1(v_vec,d_mat(:,4),dv(4),'linear',d_mat(end,4));
            d5=interp1(v_vec,d_mat(:,5),dv(5),'linear',d_mat(end,5));
            d6=interp1(v_vec,d_mat(:,6),dv(6),'linear',d_mat(end,6));
            d7=interp1(v_vec,d_mat(:,7),dv(7),'linear',d_mat(end,7));
            d8=interp1(v_vec,d_mat(:,8),dv(8),'linear',d_mat(end,8));
               
        case 'spline'
            d1=ppval(d_spline1,dv(1));
            d2=ppval(d_spline2,dv(2));
            d3=ppval(d_spline3,dv(3));
            d4=ppval(d_spline4,dv(4));
            d5=ppval(d_spline5,dv(5));
            d6=ppval(d_spline6,dv(6));
            d7=ppval(d_spline7,dv(7));
            d8=ppval(d_spline8,dv(8));
    
    end   
    D1=fct_D(d1,d2,d3,d4,d5,d6,d7,d8);
end

function D2 = dderivfunc(u_n_p_1,u_n,delta_t,type)      
    v=1/delta_t.*(u_n_p_1-u_n);
    v1=v(1);
    v2=v(2);
    v3=v(3);
    v4=v(4);
    v5=v(5);
    v6=v(6);
    v7=v(7);
    v8=v(8);
    v9=v(9);
    v10=v(10);
    v11=v(11);
    
    dv=dufunc(v);

    tmp=readmatrix('V_Fd_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    v_vec=tmp(:,1);
    d_mat=tmp(:,2:end)./v_vec;
    clear tmp;
       
    d_spline1 = spline(v_vec,d_mat(:,1));
    d_spline2 = spline(v_vec,d_mat(:,2));
    d_spline3 = spline(v_vec,d_mat(:,3));
    d_spline4 = spline(v_vec,d_mat(:,4));
    d_spline5 = spline(v_vec,d_mat(:,5));
    d_spline6 = spline(v_vec,d_mat(:,6));
    d_spline7 = spline(v_vec,d_mat(:,7));
    d_spline8 = spline(v_vec,d_mat(:,8));
    d_der1=fnder(d_spline1,1);
    d_der2=fnder(d_spline2,1);
    d_der3=fnder(d_spline3,1);
    d_der4=fnder(d_spline4,1);
    d_der5=fnder(d_spline5,1);
    d_der6=fnder(d_spline6,1);
    d_der7=fnder(d_spline7,1);
    d_der8=fnder(d_spline8,1);

    
%     global X d_grid size_grid
%     global dDdxv
%     global l_long_fl l_long_fr l_long_rl l_long_rr l_lat_fl l_lat_fr l_lat_rl l_lat_rr 
%     global l1 l2 l3 l4 l5 l6 l7 l8
%     global d_der1 d_der2 d_der3 d_der4 d_der5 d_der6 d_der7 d_der8; 
    
    switch type
        case 'linear'
            d1=interp1(v_vec,d_mat(:,1),dv(1),'linear',d_mat(end,1));
            d2=interp1(v_vec,d_mat(:,2),dv(2),'linear',d_mat(end,2));
            d3=interp1(v_vec,d_mat(:,3),dv(3),'linear',d_mat(end,3));
            d4=interp1(v_vec,d_mat(:,4),dv(4),'linear',d_mat(end,4));
            d5=interp1(v_vec,d_mat(:,5),dv(5),'linear',d_mat(end,5));
            d6=interp1(v_vec,d_mat(:,6),dv(6),'linear',d_mat(end,6));
            d7=interp1(v_vec,d_mat(:,7),dv(7),'linear',d_mat(end,7));
            d8=interp1(v_vec,d_mat(:,8),dv(8),'linear',d_mat(end,8));
            d_d1=dderivfun(d1,dv(1),1);
            d_d2=dderivfun(d2,dv(2),2);
            d_d3=dderivfun(d3,dv(3),3);
            d_d4=dderivfun(d4,dv(4),4);
            d_d5=dderivfun(d5,dv(5),5);
            d_d6=dderivfun(d6,dv(6),6);
            d_d7=dderivfun(d7,dv(7),7);
            d_d8=dderivfun(d8,dv(8),8);
            
        case 'spline'
            d_d1=ppval(d_der1,dv(1));
            d_d2=ppval(d_der2,dv(2));
            d_d3=ppval(d_der3,dv(3));
            d_d4=ppval(d_der4,dv(4));
            d_d5=ppval(d_der5,dv(5));
            d_d6=ppval(d_der6,dv(6));
            d_d7=ppval(d_der7,dv(7));
            d_d8=ppval(d_der8,dv(8));
    
    end  
     
    D2=fct_dD(d_d1,d_d2,d_d3,d_d4,d_d5,d_d6,d_d7,d_d8,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11); 
end

function dd_du = dderivfun(d,l,idx)
%     global X d_grid size_grid 

    tmp=readmatrix('V_Fd_lookup3.dat');
    tmp(tmp(:,1)==0,:)=[];
    v_vec=tmp(:,1);
    d_mat=tmp(:,2:end)./v_vec;
    clear tmp;

    
    d_grid_tmp=d_mat(:,idx);
    dnext=interp1(v_vec,d_grid_tmp,l,'next',d_grid_tmp(end));

    if d==dnext && d_grid_tmp(1)==d_grid_tmp(2)
        dd_du=0;
    elseif d==dnext && d_grid_tmp(1)~= d_grid_tmp(2) 
        dnext=d_grid_tmp(find(d_grid_tmp==d)+1);
        dprev=d_grid_tmp(find(d_grid_tmp==d)-1);
        lnext=v_vec(find(d_grid_tmp==dnext));
        lprev=v_vec(d_grid_tmp==dprev);
        dd_du=(dnext-dprev)/(lnext-lprev); 
    else
        lnext=X(d_grid_tmp==dnext);
        dd_du=(dnext-d)/(lnext-l);
    end   
end

function l=dufunc(x)
    l_long_fl=1.395;
    l_long_fr=1.395;
    l_long_rl=1.596;
    l_long_rr=1.596;
    l_lat_fl=1.6916;
    l_lat_fr=1.6916;
    l_lat_rl=1.68;
    l_lat_rr=1.68;    
    
    l1=x(1)-x(4)-l_lat_fl*x(2)-l_long_fl*x(3);
    l2=x(4)-x(5);
    l3=x(1)-x(6)+l_lat_fr*x(2)-l_long_fr*x(3);
    l4=x(6)-x(7);
    l5=x(1)-x(8)-l_lat_rl*x(2)+l_long_rl*x(3);
    l6=x(8)-x(9);
    l7=x(1)-x(10)+l_lat_rr*x(2)+l_long_rr*x(3);
    l8=x(10)-x(11);
    
    l=[l1,l2,l3,l4,l5,l6,l7,l8];
    
end

