function [F] = circular_force(pos, C, radius, M, vel)
%{
pos -> Position vector of the components of the car in the following order 
vel -> Velocity vector of the components of the car in the following order
 pos  =    [
            CG; 
            W1;
            W2;
            W3;
            W4;
            T1;
            T2;
            T3;
            T4;
            ]
C  -> Coordinate of Center of mass of the car
  C  =  [
         x;
         y;
         z
         ]

M  -> Mass of the components of the car in the following order
 M =  [
        M_cg;
        M_w1;
        M_w2;
        M_w3;
        M_w4;
        M_t1;
        M_t2;
        M_t3;
        M_t4;
        ]
For any vector in x-y plane with form a*i + b*j perpendicular in same plane
is always of form -b*i + a*j or b*i - a*j
%}
    UP = @(V)([-V(2); V(1)]);
    P_cg = pos(1:3);
    P_w1 = pos(4:6);
    P_w2 = pos(7:9);
    P_w3 = pos(10:12);
    P_w4 = pos(13:15);
    P_t1 = pos(16:18);
    P_t2 = pos(19:21);
    P_t3 = pos(22:24);
    P_t4 = pos(25:27);
    
    V_cg = vel(1:3);
    V_w1 = vel(4:6);
    V_w2 = vel(7:9);
    V_w3 = vel(10:12);
    V_w4 = vel(13:15);
    V_t1 = vel(16:18);
    V_t2 = vel(19:21);
    V_t3 = vel(22:24);
    V_t4 = vel(25:27);
    
    M_cg = M(1);
    M_w1 = M(2);
    M_w2 = M(3);
    M_w3 = M(4);
    M_w4 = M(5);
    M_t1 = M(6);
    M_t2 = M(7);
    M_t3 = M(8);
    M_t4 = M(9);

    R_cg = P_cg-C;
    R_w1 = P_w1-C;
    R_w2 = P_w2-C;
    R_w3 = P_w3-C;
    R_w4 = P_w4-C;
    R_t1 = P_t1-C;
    R_t2 = P_t2-C;
    R_t3 = P_t3-C;
    R_t4 = P_t4-C;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Check for the constraint %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eps = 1e-10;
if abs(norm(R_cg(1:2)) - radius)>eps
    error('Center of car is not on circle anymore, distance from center = %f', norm(P_cg(1:2)-C(1:2)))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Compute Direction of the Force %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CG_dir = P_cg(1:2) - C(1:2); % CG
CG_dir = CG_dir/norm(CG_dir);
w1_dir = P_w1(1:2) - C(1:2); % wheel 1
w1_dir = w1_dir / norm(w1_dir);
w2_dir = P_w2(1:2) - C(1:2); % wheel 2
w2_dir = w2_dir / norm(w2_dir);
w3_dir = P_w3(1:2) - C(1:2); % wheel 3
w3_dir = w3_dir / norm(w3_dir);
w4_dir = P_w4(1:2) - C(1:2); % wheel 4
w4_dir = w4_dir / norm(w4_dir);
t1_dir = P_t1(1:2) - C(1:2); % tire 1
t1_dir = t1_dir / norm(t1_dir);
t2_dir = P_t2(1:2) - C(1:2); % tire 2
t2_dir = t2_dir / norm(t2_dir);
t3_dir = P_t3(1:2) - C(1:2); % tire 3
t3_dir = t3_dir / norm(t3_dir);
t4_dir = P_t4(1:2) - C(1:2); % tire 4
t4_dir = t4_dir / norm(t4_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Final Force Computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vel_mag = UP(R_cg(1:2))'* V_cg(1:2) / norm(R_cg(1:2));
F_cg = (M_cg * vel_mag * vel_mag / norm(R_cg(1:2))) * CG_dir;
vel_mag = UP(R_w1(1:2))'* V_w1(1:2) / norm(R_w1(1:2));
F_w1 = (M_w1 * vel_mag* vel_mag / norm(R_w1(1:2))) * w1_dir;
vel_mag = UP(R_w2(1:2))'* V_w2(1:2) / norm(R_w2(1:2));
F_w2 = (M_w2 * (vel_mag)* vel_mag / norm(R_w2(1:2))) * w2_dir;
vel_mag = UP(R_w3(1:2))'* V_w3(1:2) / norm(R_w3(1:2));
F_w3 = (M_w3 * (vel_mag)* vel_mag / norm(R_w3(1:2))) * w3_dir;
vel_mag = UP(R_w4(1:2))'* V_w4(1:2) / norm(R_w4(1:2));
F_w4 = (M_w4 * (vel_mag)* vel_mag / norm(R_w4(1:2))) * w4_dir;
vel_mag = UP(R_t1(1:2))'* V_t1(1:2) / norm(R_t1(1:2));
F_t1 = (M_t1 * (vel_mag)* vel_mag / norm(R_t1(1:2))) * t1_dir;
vel_mag = UP(R_t2(1:2))'* V_t2(1:2) / norm(R_t2(1:2));
F_t2 = (M_t2 * (vel_mag)* vel_mag / norm(R_t2(1:2))) * t2_dir;
vel_mag = UP(R_t3(1:2))'* V_t3(1:2) / norm(R_t3(1:2));
F_t3 = (M_t3 * (vel_mag)* vel_mag / norm(R_t3(1:2))) * t3_dir;
vel_mag = UP(R_t4(1:2))'* V_t4(1:2) / norm(R_t4(1:2));
F_t4 = (M_t4 * (vel_mag)* vel_mag / norm(R_t4(1:2))) * t4_dir;

F = [F_cg;...
    F_w1;...
    F_w2;...
    F_w3;...
    F_w4;...
    F_t1;...
    F_t2;...
    F_t3;...
    F_t4];
end