L = 14;
theta = pi/4;
alpha = 1;
marker_color1 =  [019, 103, 131]/256;
marker_color2 =  [255,158,002] / 256;

% filename = ['EnergyOBCLocHamN', num2str(L), 'theta',num2str(theta, '%.4f'),'.txt'];
filename = ['EnergyOBCN', num2str(L), 'theta',num2str(theta,  '%.4f'),'alpha', num2str(alpha),'.txt'];
E_list = importdata(filename);
E0 = min(min(E_list));
E_list = E_list - E0;

% filename2 = ['BZZOBCLocHamN', num2str(L), 'theta',num2str(theta, '%.4f'),'.txt'];
filename2 = ['BZZOBCN', num2str(L), 'theta',num2str(theta,  '%.4f'),'alpha', num2str(alpha),'.txt'];
bzz_data = importdata(filename2);

% Create a logical index for positive and negative values in bzz_data
positive_idx = bzz_data > 0;
negative_idx = bzz_data < 0;

% Create markers and colors based on the logical index
markers = {'o', '+'};
colors = [marker_color1; marker_color2];

% Plot E_list with different markers and colors
hold on;
h_list = gobjects(2); % Pre-allocate array for handles
h1 = plot(-1*ones(1,sum(positive_idx(1,:))), E_list(1, positive_idx(1,:)), markers{1}, 'MarkerSize',10, 'Color', colors(1,:));
plot(1*ones(1,sum(positive_idx(2,:))), E_list(2, positive_idx(2,:)), markers{1}, 'MarkerSize',10, 'Color', colors(1,:));

h2 = plot(-1*ones(1,sum(negative_idx(1,:))),E_list(1, negative_idx(1,:)), markers{2}, 'MarkerSize',10, 'Color', colors(2,:));
plot(1*ones(1,sum(negative_idx(2,:))),E_list(2, negative_idx(2,:)), markers{2}, 'MarkerSize',10, 'Color', colors(2,:));

% % Set marker and line properties
% for h = h_list
%     set(h, 'MarkerFaceColor', get(h, 'Color'));
%     set(h, 'MarkerEdgeColor', get(h, 'Color'));
% end
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$P$','Interpreter','latex');
ylabel('$\Delta$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

grid on;
legend([h1,h2], {'# DW = even', '# DW = odd'}, 'Location', 'northwest');