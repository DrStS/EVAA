function [final_displacement] = BDF2_test(tend, h)
%clc;
%clear all;
%close all;
%format long e;
% Dynamic Backward Euler problem 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_body_fl=28e3*0.69;
k_tyre_fl=260e3;
k_body_fr=28e3*0.69;
k_tyre_fr=260e3;
k_body_rl=16e3*0.82;
k_tyre_rl=260e3;
k_body_rr=16e3*0.82;
k_tyre_rr=260e3;
l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=2*0.8458;
l_lat_fr=2*0.8458;
l_lat_rl=2*0.84;
l_lat_rr=2*0.84;
mass_Body=1936;
I_body_xx=6400;
I_body_yy=4800;
mass_wheel_fl=145/2;
mass_tyre_fl=30;
mass_wheel_fr=145/2;
mass_tyre_fr=30;
mass_wheel_rl=135/2;
mass_tyre_rl=30;
mass_wheel_rr=135/2;
mass_tyre_rr=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tend=100;
u_init=0;
du_init=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M=diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K=[k_body_fl+k_body_fr+k_body_rl+k_body_rr, k_body_fl*l_lat_fl-k_body_fr*l_lat_fr+k_body_rl*l_lat_rl-k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0  l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, -l_long_fl*l_lat_fl*k_body_fl+l_lat_fr*l_long_fr*k_body_fr+l_long_rl*l_lat_rl*k_body_rl-l_long_rr*l_lat_rr*k_body_rr, -l_lat_fl*k_body_fl, 0, l_lat_fr*k_body_fr, 0, -l_lat_rl*k_body_rl, 0, l_lat_rr*k_body_rr, 0;
   0 0 l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0 l_long_fr*k_body_fr 0 -l_long_rl*k_body_rl 0 -l_long_rr*k_body_rr 0; 
   0 0 0 k_body_fl+k_tyre_fl -k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 k_tyre_fl 0 0 0 0 0 0;
   0 0 0 0 0 k_body_fr+k_tyre_fr -k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 k_tyre_fr 0 0 0 0;
   0 0 0 0 0 0 0 k_body_rl+k_tyre_rl -k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 k_tyre_rl 0 0;
   0 0 0 0 0 0 0 0 0 k_body_rr+k_tyre_rr -k_tyre_rr;
   0 0 0 0 0 0 0 0 0 0 k_tyre_rr];
K=K+K'-diag(diag(K));
D = K *0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduced System due to Dirichlet BC (fixed to road)
Kred = K;
Kred(:,11) = [];
Kred(11,:) = [];
Kred(:,9) = [];
Kred(9,:) = [];
Kred(:,7) = [];
Kred(7,:) = [];
Kred(:,5) = [];
Kred(5,:) = [];
Mred = M;
Mred(:,11) = [];
Mred(11,:) = [];
Mred(:,9) = [];
Mred(9,:) = [];
Mred(:,7) = [];
Mred(7,:) = [];
Mred(:,5) = [];
Mred(5,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
dim_system = length(M);
t=zeros(1000+1,1);
u_sol=zeros(1000+1,11);
u_n_p_1=zeros(dim_system,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%h=1/(1000);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u_n=[u_init; zeros(dim_system - 1,1)];
u_n_m_1=u_n-h*[du_init; zeros(dim_system-1,1)];
u_n_m_2=zeros(dim_system,1);
u_n_m_3=zeros(dim_system,1);
A=((1/(h*h))*M+(1/h)*D+K);
B=((2/(h*h))*M+(1/h)*D);
C = zeros(dim_system,dim_system);
D = zeros(dim_system,dim_system);
E = zeros(dim_system,dim_system);
f_n_p_1=[1.1e3; zeros(dim_system-1,1)];
j=2;
i=h;
% first euler step
tic
if j==2 && i<tend
    rhs = (B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
    u_n_p_1 = A\rhs;
    % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;
    u_n_m_2 = u_n_m_1;
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
    i = i+h;
end
ref_val = [5.681627821952628e-07;...
     2.417538523599973e-32;...
     1.422443170437330e-12;...
     1.508286084852408e-10;...
     1.295949709390701e-12;...
     1.508286084852408e-10;...
     1.295949709390701e-12;...
     1.099930357207625e-10;...
     9.450822632980253e-13;...
     1.099930357207624e-10;...
     9.450822632980251e-13];
%disp(norm(u_n_p_1-ref_val))
% % second euler step
if j==3 && i<tend
    rhs = (B*u_n-((1/(h*h))*M)*u_n_m_1+f_n_p_1);
    u_n_p_1 = A\rhs;
    % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;
    u_n_m_3=u_n_m_2;
    u_n_m_2=u_n_m_1;
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
    i = i+h;
end

ref_val = [1.704450284873314e-06;...
    -1.695580977941673e-31;...
     7.110573430459281e-12;...
     7.529925583571522e-10;...
     9.039492654867869e-12;...
     7.529925583571522e-10;...
     9.039492654867869e-12;...
     5.490865746803850e-10;...
     6.591786028247906e-12;...
     5.490865746803851e-10;...
     6.591786028247908e-12];
%disp(norm(u_n_p_1-ref_val))


A = (9/(4*h*h))*M + (3/(2*h))*D + K;
B = (6/(h*h))*M + (2/(h))*D;
C = (-11/(2*h*h))*M + (-1/(2*h))*D;
D = (2/(h*h))*M;
E = (-1/(4*h*h))*M;
for start=i:h:tend
    rhs = B*u_n + C*u_n_m_1 + D*u_n_m_2 + E*u_n_m_3 + f_n_p_1;
    u_n_p_1 = A\rhs;
    t(j)=start;   
    u_sol(j,:)=u_n_p_1;
    u_n_m_3=u_n_m_2;
    u_n_m_2=u_n_m_1;
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
end
toc
final_displacement = u_n_p_1;
if tend==10
    ref_val = [2.354848297477015e+01;...
        -4.672535308216389e-12;...
        -1.043337042893717e-01;...
         2.369150496933984e+01;...
         2.369145008337897e+01;...
         2.369150496933973e+01;...
         2.369145008337908e+01;...
         2.337889552119290e+01;...
         2.337884792289080e+01;...
         2.337889552121592e+01;...
         2.337884792291431e+01];
    statement = sprintf('Diff for 10 sec %f',norm(u_n_p_1-ref_val));
 %   disp(statement)
elseif tend == 1
    ref_val = [2.361027257085025e-01;...
     9.086802443526854e-16;...
    -1.236265017564796e-03;...
     2.336063480730669e-01;...
     2.335138702132668e-01;...
     2.336063480730389e-01;...
     2.335138702132384e-01;...
     2.337011935281506e-01;...
     2.336951537038742e-01;...
     2.337011935281898e-01;...
     2.336951537039128e-01];
    statement = sprintf('Diff for 1 sec %f',norm(u_n_p_1-ref_val));
  %  disp(statement)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
%figure();
%plot(t,u_sol(:,1)); grid on;
%plot(t,u_sol_red(:,1)); grid on;
%legend;
%disp(u_sol(end,1:11)')
end