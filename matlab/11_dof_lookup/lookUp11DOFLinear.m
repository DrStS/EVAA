clear; clc; close all;
format long e;
%%
a = 1;
b = 2;
c = 3;

%%
global l_long_fl;
global l_long_fr;
global l_long_rl;
global l_long_rr;
global l_lat_fl;
global l_lat_fr;
global l_lat_rl;
global l_lat_rr;
l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=2*0.8458;
l_lat_fr=2*0.8458;
l_lat_rl=2*0.84;
l_lat_rr=2*0.84;
mass_Body=1936;
I_body_xx=640;
I_body_yy=4800;
mass_wheel_fl=145/2;
mass_tyre_fl=0;
mass_wheel_fr=145/2;
mass_tyre_fr=0;
mass_wheel_rl=135/2;
mass_tyre_rl=0;
mass_wheel_rr=135/2;
mass_tyre_rr=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpolation = "linear";
tend=10;
u_init=0;
du_init=0;
%%
% evaluation density
dl = 0.1;
dw = 0.1;
% calc min and max length to evaluate (just for now)
l_min = -3;
l_max = 3;
w_min = -1;
w_max = 1;
% grid size
l_d = round((l_max - l_min)/dl);
w_d = round((w_max - w_min)/dw);
% grid value allocation
global k_grid;
global X;
X = l_min:dl:l_max;
k_grid = zeros(size(X));
% k response function
k = @(l)(c*l*l + b*l + a);
% fill in grid values
for i = 0:l_d
    k_grid(i+1)= k(l_min+i*dl);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steps = 1000;
h=1/steps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate memory
t=zeros(steps+1,1);
u_sol=zeros(steps+1,11);
u_n_p_1=zeros(11,1);
u_n=[u_init; zeros(10,1)];
u_n_m_1=[u_init; zeros(10,1)]-h*[du_init; zeros(10,1)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System Mono
M=diag([mass_Body, I_body_xx, I_body_yy, mass_wheel_fl, mass_tyre_fl, mass_wheel_fr, mass_tyre_fr, mass_wheel_rl, mass_tyre_rl, mass_wheel_rr, mass_tyre_rr]);
K = get_K(u_n);
D = K *0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_pre=(1/(h*h))*M+(1/h)*D;
B=((2/(h*h))*M+(1/h)*D);
f_n_p_1=[0.1; zeros(10,1)];

%%
syms l1 l2 u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11;
l4 = u4 - u5;
l5 = u5;
l6 = u6 - u7;
l7 = u7;
l8 = u8 - u9;
l9 = u9;
l10 = u10 - u11;
l11 = u11;

%% Internal forces
syms p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 k1 k2 k3 k4 k5 k6 k7 k8 k9 k10 k11;
p1 = k4 * ( u1 - u4) + k6 * (u1 - u5) + k8 * (u1 - u8) + k10 * (u1 - u10);
% p2 = ;
% p3 = ;
p4 = k4*u4 - k5*u5;
p5 = -k4*u4+(k4+k5)*u5;
p6 = k6*u6 - k7*u7;
p7 = -k6*u6+(k6+k7)*u7;
p8 = k8*u8 - k9*u9;
p9 = -k8*u8+(k8+k9)*u9;
p10 = k10*u10 - k11*u11;
p11 = -k10*u10+(k10+k11)*u11;
%% Jacobian
syms K k_body_fl k_tyre_fl k_body_fr k_tyre_fr k_body_rl k_tyre_rl k_body_rr k_tyre_rr
% actuall derivative but sym cant be used as we get part of derivative out
% my function
k11 = diff(p1, u1);
k12 = diff(p1, u2);
k13 = diff(p1, u3);
k14 = diff(p1, u4);
k15 = diff(p1, u5);
k16 = diff(p1, u6);
k17 = diff(p1, u7);
k18 = diff(p1, u8);
k19 = diff(p1, u9);
k110 = diff(p1, u10);
k0111 = diff(p1, u11);
k21 = diff(p2, u1);
k22 = diff(p2, u2);
k23 = diff(p2, u3);
k24 = diff(p2, u4);
k25 = diff(p2, u5);
k26 = diff(p2, u6);
k27 = diff(p2, u7);
k28 = diff(p2, u8);
k29 = diff(p2, u9);
k210 = diff(p2, u10);
k211 = diff(p2, u11);
k31 = diff(p3, u1);
k32 = diff(p3, u2);
k33 = diff(p3, u3);
k34 = diff(p3, u4);
k35 = diff(p3, u5);
k36 = diff(p3, u6);
k37 = diff(p3, u7);
k38 = diff(p3, u8);
k39 = diff(p3, u9);
k310 = diff(p3, u10);
k311 = diff(p3, u11);
k41 = diff(p4, u1);
k42 = diff(p4, u2);
k43 = diff(p4, u3);
k44 = diff(p4, u4);
k45 = diff(p4, u5);
k46 = diff(p4, u6);
k47 = diff(p4, u7);
k48 = diff(p4, u8);
k49 = diff(p4, u9);
k410 = diff(p4, u10);
k411 = diff(p4, u11);
k51 = diff(p5, u1);
k52 = diff(p5, u2);
k53 = diff(p5, u3);
k54 = diff(p5, u4);
k55 = diff(p5, u5);
k56 = diff(p5, u6);
k57 = diff(p5, u7);
k58 = diff(p5, u8);
k59 = diff(p5, u9);
k510 = diff(p5, u10);
k511 = diff(p5, u11);
k61 = diff(p6, u1);
k62 = diff(p6, u2);
k63 = diff(p6, u3);
k64 = diff(p6, u4);
k65 = diff(p6, u5);
k66 = diff(p6, u6);
k67 = diff(p6, u7);
k68 = diff(p6, u8);
k69 = diff(p6, u9);
k610 = diff(p6, u10);
k611 = diff(p6, u11);
k71 = diff(p7, u1);
k72 = diff(p7, u2);
k73 = diff(p7, u3);
k74 = diff(p7, u4);
k75 = diff(p7, u5);
k76 = diff(p7, u6);
k77 = diff(p7, u7);
k78 = diff(p7, u8);
k79 = diff(p7, u9);
k710 = diff(p7, u10);
k711 = diff(p7, u11);
k81 = diff(p8, u1);
k82 = diff(p8, u2);
k83 = diff(p8, u3);
k84 = diff(p8, u4);
k85 = diff(p8, u5);
k86 = diff(p8, u6);
k87 = diff(p8, u7);
k88 = diff(p8, u8);
k89 = diff(p8, u9);
k810 = diff(p8, u10);
k811 = diff(p8, u11);
k91 = diff(p9, u1);
k92 = diff(p9, u2);
k93 = diff(p9, u3);
k94 = diff(p9, u4);
k95 = diff(p9, u5);
k96 = diff(p9, u6);
k97 = diff(p9, u7);
k98 = diff(p9, u8);
k99 = diff(p9, u9);
k910 = diff(p9, u10);
k911 = diff(p9, u11);
k101 = diff(p10, u1);
k102 = diff(p10, u2);
k103 = diff(p10, u3);
k104 = diff(p10, u4);
k105 = diff(p10, u5);
k106 = diff(p10, u6);
k107 = diff(p10, u7);
k108 = diff(p10, u8);
k109 = diff(p10, u9);
k1010 = diff(p10, u10);
k1011 = diff(p10, u11);
k111 = diff(p11, u1);
k112 = diff(p11, u2);
k113 = diff(p11, u3);
k114 = diff(p11, u4);
k115 = diff(p11, u5);
k116 = diff(p11, u6);
k117 = diff(p11, u7);
k118 = diff(p11, u8);
k119 = diff(p11, u9);
k1110 = diff(p11, u10);
k1111 = diff(p11, u11);

K_symb_1 = [k11 k12 k13 k14 k15 k16 k17 k18 k19 k110 k0111;
    k21 k22 k23 k24 k25 k26 k27 k28 k29 k210 k211;
    k31 k32 k33 k34 k35 k36 k37 k38 k39 k310 k311;
    k41 k42 k43 k44 k45 k46 k47 k48 k49 k410 k411;
    k51 k52 k53 k54 k55 k56 k57 k58 k59 k510 k511;
    k61 k62 k63 k64 k65 k66 k67 k68 k69 k610 k611;
    k71 k72 k73 k74 k75 k76 k77 k78 k79 k710 k711;
    k81 k82 k83 k84 k85 k86 k87 k88 k89 k810 k811;
    k91 k92 k93 k94 k95 k96 k97 k98 k99 k910 k911;
    k101 k102 k103 k104 k105 k106 k107 k108 k109 k1010 k1011;
    k111 k112 k113 k114 k115 k116 k117 k118 k119 k1110 k1111];

syms k_body_fl k_tyre_fl k_body_fr k_tyre_fr k_body_rl k_tyre_rl k_body_rr k_tyre_rr;
syms K


    K = [k_body_fl+k_body_fr+k_body_rl+k_body_rr, k_body_fl*l_lat_fl-k_body_fr*l_lat_fr+k_body_rl*l_lat_rl-k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0, l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, -l_long_fl*l_lat_fl*k_body_fl+l_lat_fr*l_long_fr*k_body_fr+l_long_rl*l_lat_rl*k_body_rl-l_long_rr*l_lat_rr*k_body_rr, -l_lat_fl*k_body_fl, 0, l_lat_fr*k_body_fr, 0, -l_lat_rl*k_body_rl, 0, l_lat_rr*k_body_rr, 0;
   0, 0, l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0, l_long_fr*k_body_fr, 0, -l_long_rl*k_body_rl, 0, -l_long_rr*k_body_rr, 0; 
   0, 0, 0, k_body_fl+k_tyre_fl, -k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k_body_fr+k_tyre_fr, -k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k_body_rl+k_tyre_rl, -k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k_body_rr+k_tyre_rr, -k_tyre_rr;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rr];
    K=K+K'-diag(diag(K));



f_newton = @(y_curr,y1,y2,K)( ( (1/(h*h))*M + K ) * y_curr - 2 * (1/(h*h))*M * y1 + (1/(h*h))*M * y2 - rhs);


%% Time loop
j=2;
tic
for i = h:h:tend
    K = get_K(u_n);
    A=A_pre+K;
    rhs = B*u_n-((1/(h*h))*M)*u_n_m_1 + f_n_p_1;
    u_n_p_1=A\rhs;
   % Get solution
    t(j)=i;   
    u_sol(j,:)=u_n_p_1;   
    u_n_m_1=u_n;
    u_n    =u_n_p_1;
    j=j+1;
end
toc

% Some plots
figure();
plot(t,u_sol(:,1)); grid on;
legend;


function K = get_K(x)
    global l_long_fl;
    global l_long_fr;
    global l_long_rl;
    global l_long_rr;
    global l_lat_fl;
    global l_lat_fr;
    global l_lat_rl;
    global l_lat_rr;
    upper_fl = x(4);
    lower_fl = x(5);
    upper_fr = x(6);
    lower_fr = x(7);
    upper_rl = x(8);
    lower_rl = x(9);
    upper_rr = x(10);
    lower_rr = x(11);
    
    global k_grid;
    global X;
    
    % for spline interp3(X,Y,Z,k_grid,x(2),long_fl,x(3),'spline')
    k_body_fl=interp1(X,k_grid,upper_fl - lower_fl);
    k_tyre_fl=interp1(X,k_grid,lower_fl);
    k_body_fr=interp1(X,k_grid,upper_fr-lower_fr);
    k_tyre_fr=interp1(X,k_grid,lower_fr);
    k_body_rl=interp1(X,k_grid,upper_rl-lower_rl);
    k_tyre_rl=interp1(X,k_grid,lower_rl);
    k_body_rr=interp1(X,k_grid,upper_rr-lower_rr);
    k_tyre_rr=interp1(X,k_grid,lower_rr);
    
    K = [k_body_fl+k_body_fr+k_body_rl+k_body_rr, k_body_fl*l_lat_fl-k_body_fr*l_lat_fr+k_body_rl*l_lat_rl-k_body_rr*l_lat_rr, -k_body_fl*l_long_fl-k_body_fr*l_long_fr+k_body_rl*l_long_rl+k_body_rr*l_long_rr,  -k_body_fl, 0, -k_body_fr, 0, -k_body_rl, 0, -k_body_rr, 0;
   0, l_lat_fl*l_lat_fl*k_body_fl+l_lat_fr*l_lat_fr*k_body_fr+l_lat_rl*l_lat_rl*k_body_rl+l_lat_rr*l_lat_rr*k_body_rr, -l_long_fl*l_lat_fl*k_body_fl+l_lat_fr*l_long_fr*k_body_fr+l_long_rl*l_lat_rl*k_body_rl-l_long_rr*l_lat_rr*k_body_rr, -l_lat_fl*k_body_fl, 0, l_lat_fr*k_body_fr, 0, -l_lat_rl*k_body_rl, 0, l_lat_rr*k_body_rr, 0;
   0, 0, l_long_fl*l_long_fl*k_body_fl+l_long_fr*l_long_fr*k_body_fr+l_long_rl*l_long_rl*k_body_rl+l_long_rr*l_long_rr*k_body_rr, l_long_fl*k_body_fl, 0, l_long_fr*k_body_fr, 0, -l_long_rl*k_body_rl, 0, -l_long_rr*k_body_rr, 0; 
   0, 0, 0, k_body_fl+k_tyre_fl, -k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, k_tyre_fl, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, k_body_fr+k_tyre_fr, -k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, k_tyre_fr, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, k_body_rl+k_tyre_rl, -k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rl, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, k_body_rr+k_tyre_rr, -k_tyre_rr;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, k_tyre_rr];
    K=K+K'-diag(diag(K));
end