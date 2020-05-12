% test stuff
x1 = rand(1,3);
y1 = rand(1,3);
x2 = rand(1,3);
y2 = rand(1,3);

close all;

% modify values here

% plot line size
line_size_dataset_1 = 2;
line_size_dataset_2 = 2;

% title entries
title_name = 'Title Name';
title_fontsize = 20;

% labels entries
xlabel_name = 'x label $[sec]$';
ylabel_name = 'y label $[fancy]$';
label_fontsize = 16;
x_lim_down = 0;
x_lim_up = 1;
y_lim_down = 0;
y_lim_up = 1;

% legend entries
legend_name_1 = '$entry_1$';
legend_name_2 = '$entry_2$';
legend_size = 17;

% where to use stuff
figure(1); hold on; grid on;

plot(x1, y1, 'LineWidth', line_size_dataset_1);
plot(x2, y2, 'LineWidth', line_size_dataset_2);

title(title_name, 'Interpreter', 'latex', 'Fontsize', title_fontsize); 

xlim([x_lim_down, x_lim_up]);
ylim([y_lim_down, y_lim_up]);

xlabel(xlabel_name, 'Interpreter', 'latex', 'Fontsize', label_fontsize);
ylabel(ylabel_name, 'Interpreter', 'latex', 'Fontsize', label_fontsize);

leg1 = legend('$entry_1$','$entry_2$');
set(leg1,'Interpreter','latex');
set(leg1,'FontSize', legend_size);
leg1.Location = 'best';

