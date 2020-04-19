clear all
clc
ALE = dlmread('ALE_result.dat');
MBD = dlmread('MBD_result.dat');
ALE_init = dlmread('initial_car_pos_vec.dat');
idx = 3:3:27;
ALE_z = ALE_init(idx);
n = size(ALE);
n = n(1);
summand = repmat(ALE_z,n,1);
ALE(:,idx) = ALE(:, idx) + summand;
mbd_idx = [34, 36, 35, 46, 48, 47, 58, 60, 59,43, 45, 44,55, 57, 56, 37, 39, 38, 49, 51, 50,40, 42, 41,52, 54, 53];
mbd_idx = mbd_idx + ones(size(mbd_idx));
MBD_mod = MBD(:, mbd_idx);
diff = MBD_mod - ALE;
nrms = vecnorm(diff,2,2);
mbd_nrms = vecnorm(MBD_mod,2,2);
rel_nrms = nrms./mbd_nrms;
t = 2:1:n;
figure(1)
hold on
plot(t , rel_nrms(2:n,:))
title("Relative Norm - ALE vs MBD")
xlabel("Time steps")
ylabel("${\|x_{ALE}-x_{MBD}\|}$ / ${\|x_{MBD}\|}$",'Interpreter','latex','FontWeight','bold','FontSize',22)
hold off

xy_idx = [34, 36, 46, 48, 58, 60, 43, 45,55, 57, 37, 39, 49, 51, 40, 42, 52, 54];
xy_idx = xy_idx + ones(size(xy_idx));
MBD_mod_xy = MBD(:, xy_idx);
xy_idx = [1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23,25,26];
ALE_xy = ALE(:, xy_idx);
diff = MBD_mod_xy - ALE_xy;
nrms = vecnorm(diff,2,2);
mbd_nrms = vecnorm(MBD_mod_xy,2,2);
rel_nrms = nrms./mbd_nrms;

t = 2:1:n;
figure(2)
hold on
plot(t , rel_nrms(2:n,:))
title("Relative Norm(X, Y) - ALE vs MBD")
xlabel("Time steps")
ylabel("${\|x_{ALE}-x_{MBD}\|}$ / ${\|x_{MBD}\|}$",'Interpreter','latex','FontWeight','bold','FontSize',22)
hold off


%%%%%%%%%%%%% z norm
z_idx = [ 35 47 59, 44, 56, 38,50, 41 53];
z_idx = z_idx + ones(size(z_idx));
MBD_mod_z = MBD(:, z_idx);
z_idx = [3,6,9,12,15,18,21,24,27];
ALE_z = ALE(:, z_idx);
diff = MBD_mod_z - ALE_z;

nrms = vecnorm(diff,2,2);
mbd_nrms = vecnorm(MBD_mod_z,2,2);
rel_nrms = nrms./mbd_nrms;
t = 2:1:n;
figure(3)
hold on
plot(t , rel_nrms(2:n,:))
title("Relative Norm(Z) - ALE vs MBD")
xlabel("Time steps")
ylabel("${\|x_{ALE}-x_{MBD}\|}$ / ${\|x_{MBD}\|}$",'Interpreter','latex','FontWeight','bold','FontSize',22)
hold off