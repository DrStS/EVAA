clc;
clear all;
close all;
format long e;
h=1/1000;
t = zeros(100,1);
j=1;
diff = zeros(100,1);
for tend = 1:1:100
    disp_euler = BEoneFile03(tend,h);
    disp_bdf = BDF2_test(tend,h);
    diff(j)= norm(disp_bdf-disp_euler);
    t(j) = tend;
    j = j+1;
end
figure();
plot(t,diff); grid on;
legend;
%statement = sprintf('Diff for %f sec %f',tend,norm(disp_bdf-disp_euler));
%disp(statement)