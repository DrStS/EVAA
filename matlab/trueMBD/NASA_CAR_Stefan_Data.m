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
penalty = 1e5;
k_body_rot_fl = penalty;  % this is a penalty value->try it as high as possible and see how the convergence evolves
k_body_rot_fr = penalty;
k_body_rot_rl = penalty;
k_body_rot_rr = penalty;
k_tyre_rot_fl = penalty;
k_tyre_rot_fr = penalty;
k_tyre_rot_rl = penalty;
k_tyre_rot_rr = penalty;  % 
c_body_fl=@(v)0.2*k_body_fl*v;    % allow nonlinear damping 
c_tyre_fl=@(v)0.2*k_tyre_fl*v;
c_body_fr=@(v)0.2*k_body_fr*v;
c_tyre_fr=@(v)0.2*k_tyre_fr*v;
c_body_rl=@(v)0.2*k_body_rl*v;
c_tyre_rl=@(v)0.2*k_tyre_rl*v;
c_body_rr=@(v)0.2*k_body_rr*v;
c_tyre_rr=@(v)0.2*k_tyre_rr*v;
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

% simulation specifications 
num_iter = 10000;     % LENGTH OF THE SIMULATION
delta_t = 1e-3;       % PLAY AROUND
tol = 1e-7;           % !!!! PLAY AROUND !!!!
max_iter = 1000000;     % for Broyden
 
% Dimensions of the main car body (the center of rotation is at the origin)
r1 = [-l_long_rr; 0  ;  l_lat_rr];
r2 = [-l_long_rl; 0  ; -l_lat_rl];
r3 = [l_long_fl ; 0  ; -l_lat_fl];
r4 = [l_long_fr ; 0  ;  l_lat_fr];

mass = mass_Body;
Ic = diag([I_body_xx, I_body_yy, I_body_zz]);                           % moment of intertia of the car

% wheel parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_wheel = [mass_wheel_rr, mass_wheel_rl, mass_wheel_fl, mass_wheel_fr];      

% tyre parameters (provided as vectors [right-back, left-back, left-front, right_front])
mass_tyre = [mass_tyre_rr, mass_tyre_rl, mass_tyre_fl, mass_tyre_fr];                   % do not put to zero to avoid singularities, Dirichlet condition are enforced via a force vector

upper_spring_length = [ 0.7889; 0.7889; 0.7069; 0.7069];
lower_spring_length = [ 0.3; 0.3; 0.3; 0.3];

initial_upper_spring_length = [0.4470; 0.4470; 0.4470; 0.4470];                 % PLAY AROUND 
initial_lower_spring_length = [0.2769; 0.2769; 0.2769; 0.2769];

upper_spring_stiffness = [k_body_rr; k_body_rl; k_body_fl; k_body_fr];
lower_spring_stiffness = [k_tyre_rr; k_tyre_rl; k_tyre_fl; k_tyre_fr];

upper_spring_damping = {c_body_rr, c_body_rl, c_body_fl, c_body_fr};
lower_spring_damping = {c_tyre_rr, c_tyre_rl, c_tyre_fl, c_tyre_fr};

upper_rotational_stiffness = [k_body_rot_rr; k_body_rot_rl; k_body_rot_fl; k_body_rot_fr];
lower_rotational_stiffness = [k_tyre_rot_rr; k_tyre_rot_rl; k_tyre_rot_fl; k_tyre_rot_fr];

generate_from_trajectoryXYZ = true;

% allocalte memory
F_tyre1 = zeros(3,num_iter+1);
F_tyre2 = zeros(3,num_iter+1);
F_tyre3 = zeros(3,num_iter+1);
F_tyre4 = zeros(3,num_iter+1);


if generate_from_trajectoryXYZ
    % create trajectory
%     v_init = 3;    
%     radius = 5;
%     trajectoryXYZ = generate_circular_trajectory(radius, initial_position, initial_direction, delta_t, v_init, num_iter+1);
    
%    v_init = 3;    
%    amplitude = 2;
%    period = 20;

%    initial_position = [10; 40; 20];
%    initial_direction = [1;-1];
%    trajectoryXYZ = generate_drunk_car(amplitude, period, initial_position, initial_direction, delta_t, v_init, num_iter + 1);

%     acceleration = 2;
%     v_init = 0.1;
%     initial_position = [10; 40; 20];
%     initial_direction = [1;0];
%     trajectoryXYZ = generate_accelerated_car(acceleration, initial_position, initial_direction, delta_t, v_init, num_iter + 1);

%     ramp_length = 60;
%     pre_path = 10;
%     layout = @(X)0.001 * X^2;
%     initial_position = [10; 40; 20];
%     initial_direction = [-1;0];
%     jump = -2;
%     v_init = 10;
%     trajectoryXYZ = generate_ramp(pre_path, ramp_length, layout, initial_position, initial_direction, jump, delta_t, v_init, num_iter + 1);
    
   [trajectoryXYZ, num_iter] = generate_fancy_road(delta_t);

    % get the trajectories of the tyres
    [trajectoryXYZ_t1, trajectoryXYZ_t2, trajectoryXYZ_t3, trajectoryXYZ_t4, theta, vc, w] ...
        = compute_all_tyre_positions([-l_long_rr; l_lat_rr],...
                                     [-l_long_rl; -l_lat_rl],...
                                     [l_long_fl ; -l_lat_fl], ...
                                     [l_long_fr ; l_lat_fr], trajectoryXYZ, delta_t);
                                 
                                 
    % get Z components of the tyre positions
    trajectoryXYZ_t1(2, :) = calculate_Zpositions_tyres(trajectoryXYZ, r1(1), -initial_upper_spring_length(1) -initial_lower_spring_length(1));
    trajectoryXYZ_t2(2, :) = calculate_Zpositions_tyres(trajectoryXYZ, r2(1), -initial_upper_spring_length(2) -initial_lower_spring_length(2));
    trajectoryXYZ_t3(2, :) = calculate_Zpositions_tyres(trajectoryXYZ, r3(1), -initial_upper_spring_length(3) -initial_lower_spring_length(3));
    trajectoryXYZ_t4(2, :) = calculate_Zpositions_tyres(trajectoryXYZ, r4(1), -initial_upper_spring_length(4) -initial_lower_spring_length(4));
   
   % get the forces in the tyres
   [F_tyre1, ~] = get_forces_arbitrary_path(trajectoryXYZ_t1, theta, mass_tyre(1), Ic(3,3), delta_t);
   [F_tyre2, ~] = get_forces_arbitrary_path(trajectoryXYZ_t2, theta, mass_tyre(2), Ic(3,3), delta_t);
   [F_tyre3, ~] = get_forces_arbitrary_path(trajectoryXYZ_t3, theta, mass_tyre(3), Ic(3,3), delta_t);
   [F_tyre4, ~] = get_forces_arbitrary_path(trajectoryXYZ_t4, theta, mass_tyre(4), Ic(3,3), delta_t);

   % initial velocities
    vc = vc(:,1);       % car body

    [vt1, vt2, vt3, vt4] = get_all_path_velocities(trajectoryXYZ_t1, trajectoryXYZ_t2, trajectoryXYZ_t3, trajectoryXYZ_t4, delta_t);
    [vw1, vw2, vw3, vw4] = get_all_path_velocities(trajectoryXYZ_t1, trajectoryXYZ_t2, trajectoryXYZ_t3, trajectoryXYZ_t4, delta_t);

    % initial angular velocities
    wc = [0; w(1); 0];
    initial_orientation = calculate_quaternion_from_axis_angle([0;-1;0], theta(1));                             
    if isnan(initial_orientation)
        initial_orientation = [0; 0; 0; 1];
    end
    initial_position = trajectoryXYZ(:,1);  

else % OLD VERSION
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
    wc = [0; 0; 0];

    initial_orientation = [0; -1; 0; 1];                                     % PLAY AROUND WITH IT initial orientation of the car body as quaternion
    initial_position = [5;0;0];                                              % of the center of mass
end

visualize = true;

% force parameters
g = 9.81;               % there is no gravity in outer space! 

FC = [0; -mass*g; 0];    % external forces in y_direction

FT1 = [0; -mass_tyre(1)*g; 0];
FT2 = [0; -mass_tyre(2)*g; 0];
FT3 = [0; -mass_tyre(3)*g; 0];
FT4 = [0; -mass_tyre(4)*g; 0];

FW1 = [0; -mass_wheel(1)*g; 0];
FW2 = [0; -mass_wheel(2)*g; 0];
FW3 = [0; -mass_wheel(3)*g; 0];
FW4 = [0; -mass_wheel(4)*g; 0];


% road forces in y direction (as functions of time)
%fancy road forces (do not work right now)
% FR1 = @(t, i, y, pcc, vt, vb, F) flying_car_road_forces(y, vt, mass_tyre(1), F, F_tyre1(i), trajectoryXYZ_pt1(2,i), delta_t);        
% FR2 = @(t, i, y, pcc, vt, vb, F) flying_car_road_forces(y, vt, mass_tyre(2), F, F_tyre2(i), trajectoryXYZ_pt2(2,i), delta_t);
% FR3 = @(t, i, y, pcc, vt, vb, F) flying_car_road_forces(y, vt, mass_tyre(3), F, F_tyre3(i), trajectoryXYZ_pt3(2,i), delta_t);
% FR4 = @(t, i, y, pcc, vt, vb, F) flying_car_road_forces(y, vt, mass_tyre(4), F, F_tyre4(i), trajectoryXYZ_pt4(2,i), delta_t);

%the car is flying away 
% FR1 = @(t, y, pcc, vt, vb, F) F;        
% FR2 = @(t, i, y, pcc, vt, vb, F) F;
% FR3 = @(t, i, y, pcc, vt, vb, F) F;
% FR4 = @(t, i, y, pcc, vt, vb, F) F;

%fix the car on the ground (initial velocities very small) -> MOST SIMULATIONS WITH THIS 
% FR1 = @(t, i, y, pcc, vt, vb, F) zeros(3,1);
% FR2 = @(t, i, y, pcc, vt, vb, F) zeros(3,1);
% FR3 = @(t, i, y, pcc, vt, vb, F) zeros(3,1);
% FR4 = @(t, i, y, pcc, vt, vb, F) zeros(3,1);

% FR1 = @(t, i, y, pcc, vt, vc, F) Circular_path(vt, mass_tyre(1), y); % If you don't use circular path, uncomment line 10 in main_nasa_car
% FR2 = @(t, i, y, pcc, vt, vc, F) Circular_path(vt, mass_tyre(2), y);
% FR3 = @(t, i, y, pcc, vt, vc, F) Circular_path(vt, mass_tyre(3), y);
% FR4 = @(t, i, y, pcc, vt, vc, F) Circular_path(vt, mass_tyre(4), y);

if generate_from_trajectoryXYZ
    % arbitrary trajectory
    FR1 = @(t, i, y, pcc, vt, vb, F) ...
        flying_car_road_forces(mass_tyre(1), F, 1, F_tyre1, trajectoryXYZ_t1, vt1, i, delta_t);        
    FR2 = @(t, i, y, pcc, vt, vb, F) ...
        flying_car_road_forces(mass_tyre(2), F, 2, F_tyre2, trajectoryXYZ_t2, vt2, i, delta_t);
    FR3 = @(t, i, y, pcc, vt, vb, F) ...
        flying_car_road_forces(mass_tyre(3), F, 3, F_tyre3, trajectoryXYZ_t3, vt3, i, delta_t);
    FR4 = @(t, i, y, pcc, vt, vb, F) ...
        flying_car_road_forces(mass_tyre(4), F, 4, F_tyre4, trajectoryXYZ_t4, vt4, i, delta_t);
end

%PLAY AROUND; expect RK4 and BCN to work fine
% Explicit solvers 
% solver = @(f, t, x) explicit_solver(f, t, x);
solver = @(f, t, x) Runge_Kutta_4(f, t, x);

% Implicit solvers
% solver = @(f, t, x) Broyden_Euler(f, t, x, tol, max_iter);
% solver = @(f, t, x) Broyden_Crank_Nicolson(f, t, x, tol, max_iter);
% solver = @(f, t, x) Broyden_PDF2(f, t, x, tol, max_iter);

%% solving
[t,y, y_sol, metrics] =  main_nasa_car(r1, r2, r3, r4, mass, mass_wheel, mass_tyre, Ic, initial_orientation, initial_position, ... 
                lower_spring_length, upper_spring_length, initial_lower_spring_length,...
                trajectoryXYZ_t1,trajectoryXYZ_t2, trajectoryXYZ_t3, trajectoryXYZ_t4,  ...
                lower_spring_stiffness, upper_spring_stiffness, lower_spring_damping, upper_spring_damping, lower_rotational_stiffness, upper_rotational_stiffness, ...
                vc, vw1(:,1), vw2(:,1), vw3(:,1), vw4(:,1), vt1(:,1), vt2(:,1), vt3(:,1), vt4(:,1), wc, FC, FW1, FW2, FW3, FW4, FT1, FT2, FT3, FT4, FR1, FR2, FR3, FR4, ...
                num_iter, delta_t, solver);

%% visulalization
vis_fig = figure(1);
err_fig = figure(2);
iter_fig = figure(3);
cond_fig = figure(4);
center_fig = figure(5);

posstr=['Final y-position of the center of mass: ', num2str(y(end,6))];
disp(posstr);

set(0,'CurrentFigure',center_fig)
plot(t,y(:,6));
title("Evolution of the y-coordinate of the center of mass of the car")

if (size(metrics, 1)>1)
    set(0,'CurrentFigure',iter_fig)
    semilogy(t, metrics(:,1));
    dispstring = ["Number of Broyden iterations required to reach the tolerance ", tol];
    title(dispstring)

    set(0,'CurrentFigure',cond_fig)
    semilogy(t, metrics(:,2));
    dispstring = "Condition number of the approximated Jacobian";
    title(dispstring)
end

error_plotter3D(err_fig, solver, t, num_iter, y, y_sol, mass, mass_wheel, mass_tyre, Ic, g,...
                upper_spring_length, lower_spring_length, ...
                upper_spring_stiffness, lower_spring_stiffness, ...
                upper_rotational_stiffness, lower_rotational_stiffness);

vel_norms = zeros(1, length(t));
for i=1:length(t)
    vel_norms(i) = norm(y_sol(i, 4:6));
end

if visualize
    visualizer3D(vis_fig, y, trajectoryXYZ_t1, trajectoryXYZ_t2, delta_t, vel_norms);
end



