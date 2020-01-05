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
mass_tyre_fl=30;
mass_wheel_fr=145/2;
mass_tyre_fr=30;
mass_wheel_rl=135/2;
mass_tyre_rl=30;
mass_wheel_rr=135/2;
mass_tyre_rr=30;

% Dimensions of the main car body (the center of rotation is at the origin)
r1 = [-l_long_rr; 0  ;  l_lat_rr];
r2 = [-l_long_rl; 0  ; -l_lat_rl];
r3 = [l_long_fl ; 0  ; -l_lat_fl];
r4 = [l_long_fr ; 0  ;  l_lat_fr];

mass = mass_Body;
%center_of_mass = [length/2; 0; -width/2];                              %%NOT WORKING!!! use if rotation and mass centers are not at the same positions
Ic = diag([I_body_xx, 0, I_body_yy]);                                   % moment of intertia of the car
initial_orientation = [0; 0; 0; 1];                                     % initial orientation of the car body as quaternion

% wheel parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_wheel = [mass_wheel_rr, mass_wheel_rl, mass_wheel_fl, mass_wheel_fr];      

% tyre parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_tyre = [mass_tyre_rr, mass_tyre_rl, mass_tyre_fl, mass_tyre_fr];                   % do not put to zero to avoid singularities, Dirichlet condition are enforced via a force vector
% use reduced == true to run Stefans reduced 7DOF
reduced = false;

lower_spring_length = [0.2; 0.2; 0.2; 0.2];
upper_spring_length = [0.2; 0.2; 0.2; 0.2];

initial_lower_spring_length = [0.1; 0.1; 0.1; 0.1];
initial_upper_spring_length = [0.2; 0.2; 0.2; 0.2];

lower_spring_stiffness = [k_body_rr; k_body_rl; k_body_fl; k_body_fr];
upper_spring_stiffness = [k_tyre_rr; k_tyre_rl; k_tyre_fl; k_tyre_fr];

% initial velocities (only y-components)
vc = 0;               % car body

vw = [0; 0; 0; 0];    % wheel 

vt = [0; 0; 0; 0];    % tyres


% initial angular velocities
wc = zeros(3,1);

% force parameters
g = 0;               % there is no gravity in outer space!

%FC = -mass*g;    % external forces in y_direction
FC = 1.1e3; 

FT1 = -mass_tyre(1)*g;
FT2 = -mass_tyre(2)*g;
FT3 = -mass_tyre(3)*g;
FT4 = -mass_tyre(4)*g;

FW1 = -mass_wheel(1)*g;
FW2 = -mass_wheel(2)*g;
FW3 = -mass_wheel(3)*g;
FW4 = -mass_wheel(4)*g;

% road forces in y direction (as functions of time)
% use this to fix the car on the ground
FR1 = @(t,F) -F;        
FR2 = @(t,F) -F;
FR3 = @(t,F) -F;
FR4 = @(t,F) -F;
%the car is flying away 
%FR1 = @(t,F) 0;        
%FR2 = @(t,F) 0;
%FR3 = @(t,F) 0;
%FR4 = @(t,F) 0;

% simulation specifications
num_iter = 1e3;
delta_t = 1e-3;
tol = 1e-10;
max_iter = 200;

% Explicit solvers
% solver = @(f, t, x) explicit_solver(f, t, x);
% solver = @(f, t, x) ode45(f, t, x);
% solver = @(f, t, x) ode113(f, t, x);
% solver = @(f, t, x) Runge_Kutta_4(f, t, x);

% Implicit solvers
 solver = @(f, t, x) Broyden_Euler(f, t, x, tol, max_iter);
% solver = @(f, t, x) Broyden_Crank_Nicolson(f, t, x, tol, max_iter);
% solver = @(f, t, x) Broyden_PDF2(f, t, x, tol, max_iter);

%% solving
[t,y, y_all] = main_nasa_car(r1, r2, r3, r4, mass, mass_wheel, mass_tyre, Ic, initial_orientation, reduced,... 
        lower_spring_length, upper_spring_length, initial_lower_spring_length, initial_upper_spring_length, lower_spring_stiffness, upper_spring_stiffness,...
        vc, vw, vt, wc, FC, FW1, FW2, FW3, FW4, FT1, FT2, FT3, FT4, FR1, FR2, FR3, FR4,...
        num_iter, delta_t, solver);


   
%% error measure
error = zeros(size(t,2), 1);
                      
for i = 1 : num_iter

    %% get components
    pcc = y(i,5);

    pc = y(i,6:9);

    pw = y(i,10);

    pt = y(i,14:17);
    
    upper_length = pc - pw;
    
    lower_length = pw - pt;
    
    wc = y_all(i,1:3);
    
    vc = y_all(i,4);

    vw = y_all(i,5:8);

    vt = y_all(i,9:12);
    
    %compute the energies in the system
    potential_energy = g * (mass * pcc + sum(mass_wheel .* pw) + sum(mass_tyre .* pt));
    kinetic_energy = 0.5 * (mass * vc^2 + sum(mass_wheel .* vw.^2) + sum(mass_tyre .* vt.^2));
    rotational_energy = 0.5 * wc * Ic * wc';
    spring_potential = 0.5 * (sum(upper_spring_stiffness .* (upper_spring_length - upper_length').^2) ...
        + sum(lower_spring_stiffness .* (lower_spring_length - lower_length').^2));
    
    %compute the total energy
    error(i) = potential_energy + kinetic_energy + rotational_energy + spring_potential;
    
    % compute the relative error
    if i > 1
        error(i) = (error(i) - error(1)) / error(1);
    end
end
    error(1) = 0;

plot(t(2:end-1), error(2:end-1)); grid on;
legend('energy variation in the 3D system');
title('Error graph - reduced system, CN, 300 sec, dt=1e-3');
xlabel('Time [sec]');
ylabel('Energy');
