clear; close all; clc;

%% initial conditions
k_body_fl=28e3*0.69;
k_tyre_fl=260e3;
k_body_fr=28e3*0.69;
k_tyre_fr=260e3;
k_body_rl=16e3*0.82;
k_tyre_rl=260e3;
k_body_rr=16e3*0.82;
k_tyre_rr=260e3;
k_body_rot_fl = 1e4;
k_body_rot_fr = 1e4;
k_body_rot_rl = 1e4;
k_body_rot_rr = 1e4;
k_tyre_rot_fl = 1e4;
k_tyre_rot_fr = 1e4;
k_tyre_rot_rl = 1e4;
k_tyre_rot_rr = 1e4;
c_body_fl=@(v)10*v;    %allow nonlinear damping
c_tyre_fl=@(v)10*v;
c_body_fr=@(v)10*v;
c_tyre_fr=@(v)10*v;
c_body_rl=@(v)10*v;
c_tyre_rl=@(v)10*v;
c_body_rr=@(v)10*v;
c_tyre_rr=@(v)10*v;
l_long_fl=1.395;
l_long_fr=1.395;
l_long_rl=1.596;
l_long_rr=1.596;
l_lat_fl=2*0.84;
l_lat_fr=2*0.84;
l_lat_rl=2*0.84;
l_lat_rr=2*0.84;
mass_Body=1936;
I_body_xx=6400;
I_body_yy=4800;
I_body_zz=5000;
mass_wheel_fl=145/2;
mass_tyre_fl=30;
mass_wheel_fr=145/2;
mass_tyre_fr=30;
mass_wheel_rl=135/2;
mass_tyre_rl=30;
mass_wheel_rr=135/2;
mass_tyre_rr=30;

visualize = false;

% Dimensions of the main car body (the center of rotation is at the origin)
r1 = [-l_long_rr; 0  ;  l_lat_rr];
r2 = [-l_long_rl; 0  ; -l_lat_rl];
r3 = [l_long_fl ; 0  ; -l_lat_fl];
r4 = [l_long_fr ; 0  ;  l_lat_fr];

mass = mass_Body;
Ic = diag([I_body_xx, I_body_zz, I_body_yy]);                           % moment of intertia of the car
initial_orientation = [0; 0; 0; 1];                                     % DON'T CHANGE!! initial orientation of the car body as quaternion

% wheel parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_wheel = [mass_wheel_rr, mass_wheel_rl, mass_wheel_fl, mass_wheel_fr];      

% tyre parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_tyre = [mass_tyre_rr, mass_tyre_rl, mass_tyre_fl, mass_tyre_fr];                   % do not put to zero to avoid singularities, Dirichlet condition are enforced via a force vector


upper_spring_length = [0.2; 0.2; 0.2; 0.2];
lower_spring_length = [0.2; 0.2; 0.2; 0.2];

initial_upper_spring_length = [0.2; 0.2; 0.2; 0.2];
initial_lower_spring_length = [0.16; 0.16; 0.16; 0.16];

upper_spring_stiffness = [k_body_rr; k_body_rl; k_body_fl; k_body_fr];
lower_spring_stiffness = [k_tyre_rr; k_tyre_rl; k_tyre_fl; k_tyre_fr];

upper_spring_damping = {c_body_rr, c_body_rl, c_body_fl, c_body_fr};
lower_spring_damping = {c_tyre_rr, c_tyre_rl, c_tyre_fl, c_tyre_fr};

upper_rotational_stiffness = [k_body_rot_rr; k_body_rot_rl; k_body_rot_fl; k_body_rot_fr];
lower_rotational_stiffness = [k_tyre_rot_rr; k_tyre_rot_rl; k_tyre_rot_fl; k_tyre_rot_fr];

% initial velocities 
vc = [0; 0; 0];       % car body

vw1 = [0; 0; 0];    % wheel 
vw2 = [0; 0; 0];    
vw3 = [0; 0; 0];    
vw4 = [0; 0; 0];    

vt1 = [0; 0; 0];    % tyres 
vt2 = [0; 0; 0];    
vt3 = [0; 0; 0];    
vt4 = [0; 0; 0];    


% initial angular velocities
wc = zeros(3,1);

% force parameters
g = 0;               % there is no gravity in outer space!

FC = -mass*g;    % external forces in y_direction
%FC = [0; 1.1e3; 0]; 

FT1 = -mass_tyre(1)*g;
FT2 = -mass_tyre(2)*g;
FT3 = -mass_tyre(3)*g;
FT4 = -mass_tyre(4)*g;

FW1 = -mass_wheel(1)*g;
FW2 = -mass_wheel(2)*g;
FW3 = -mass_wheel(3)*g;
FW4 = -mass_wheel(4)*g;

% simulation specifications
num_iter = 1e5;
delta_t = 1e-3;
tol = 1e-10;
max_iter = 10000;

 
% road forces in y direction (as functions of time)
% use this to fix the car on the ground
d = -0.5;
FR1 = @(t, y, v, m, F) flying_car_road_forces(y, v, m, F, d, delta_t);        
FR2 = @(t, y, v, m, F) flying_car_road_forces(y, v, m, F, d, delta_t);
FR3 = @(t, y, v, m, F) flying_car_road_forces(y, v, m, F, d, delta_t);
FR4 = @(t, y, v, m, F) flying_car_road_forces(y, v, m, F, d, delta_t);

%the car is flying away 
FR1 = @(t, y, v, m, F) 0;        
FR2 = @(t, y, v, m, F) (t<0)*300;
FR3 = @(t, y, v, m, F) (t<0)*300;
FR4 = @(t, y, v, m, F) 0;


% Explicit solvers
% solver = @(f, t, x) explicit_solver(f, t, x);
% solver = @(f, t, x) Runge_Kutta_4(f, t, x);

% Implicit solvers
% solver = @(f, t, x) Broyden_Euler(f, t, x, tol, max_iter);
 solver = @(f, t, x) Broyden_Crank_Nicolson(f, t, x, tol, max_iter);
% solver = @(f, t, x) Broyden_PDF2(f, t, x, tol, max_iter);
%% solving
[t,y, y_sol] =  main_nasa_car(r1, r2, r3, r4, mass, mass_wheel, mass_tyre, Ic, initial_orientation,... 
                lower_spring_length, upper_spring_length, initial_lower_spring_length, initial_upper_spring_length, ...
                lower_spring_stiffness, upper_spring_stiffness, lower_spring_damping, upper_spring_damping, lower_rotational_stiffness, upper_rotational_stiffness, ...
                vc, vw1, vw2, vw3, vw4, vt1, vt2, vt3, vt4, wc, FC, FW1, FW2, FW3, FW4, FT1, FT2, FT3, FT4, FR1, FR2, FR3, FR4,...
                num_iter, delta_t, solver);

%% visulalization
figure()
plot(t,y(:,6));
posstr=['Final y-position of the center of mass: ', num2str(y(end,6))];
disp(posstr);
title("Evolution of the y-coordinate of the center of mass of the car")
if visualize
    visualizer3D(y, delta_t);
end

error_plotter3D(solver, t, num_iter, y, y_sol, mass, mass_wheel, mass_tyre, Ic, g,...
                upper_spring_length, lower_spring_length, ...
                upper_spring_stiffness, lower_spring_stiffness, ...
                upper_rotational_stiffness, lower_rotational_stiffness);


