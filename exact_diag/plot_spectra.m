L = 24;
theta = 0.6;
alpha = 1;
marker_color1 =  [019, 103, 131]/256;
marker_color2 =  [255,158,002] / 256;


% filename = ['EnergyLocHamN', num2str(L), 'theta',num2str(theta, '%.4f'),'.txt'];
filename = ['EnergyN', num2str(L), 'theta',num2str(theta, '%.4f'),'alpha',num2str(alpha), '.txt'];
data = importdata(filename);
E_list = data;
E0 = min(min(E_list));
gap1 = min(E_list(2,:))-E0;
E_list = E_list - E0;
E_list = E_list /gap1;

E_list = [E_list(L/2+1:end,:); E_list(1:L/2,:)];



E_list = E_list(:);
x_list = repmat(-L/2:L/2-1, 1, size(data, 2));


% filename2 = ['UXLocHamN', num2str(L), 'theta',num2str(theta,  '%.4f'), '.txt'];
filename2 = ['UXN', num2str(L), 'theta',num2str(theta,  '%.4f'), 'alpha', num2str(alpha),'.txt'];
ux_data = importdata(filename2);

ux_data = [ux_data(L/2+1:end,:); ux_data(1:L/2,:)];
ux_data = ux_data(:);
% Create a logical index for positive and negative values in ux_data
positive_idx = ux_data > 0;
negative_idx = ux_data < 0;


h2 = plot(x_list(negative_idx),E_list(negative_idx),'o', 'MarkerSize',10);hold on;
% set(h2, 'MarkerFaceColor', marker_color2);
set(h2, 'MarkerEdgeColor', marker_color2);

h1 = plot(x_list(positive_idx),E_list(positive_idx),'+', 'MarkerSize',10); 
% set(h1, 'MarkerFaceColor', marker_color1);
set(h1, 'MarkerEdgeColor', marker_color1);



set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$S$','Interpreter','latex');
ylabel('$\Delta$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);
grid on;

legend([h1,h2], {'$U(X)=1$', '$U(X)=-1$'},'Interpreter','latex', 'Location', 'northwest');

xlim([-inf,0]);
ax = gca;
% Get the position of the axis in normalized units (0 to 1)
pos = ax.Position;
xL=xlim;
yL=ylim;
% Create the text object
T = text(0.99*xL(2), 0.99*yL(2), ['L = ', num2str(L)], 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

set(T, 'fontsize', 24);

% Set the figure size and position
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8 6]);

% Set the filename and path
filename = ['../figure/EnergySpectraL', num2str(L), 'theta', num2str(theta, '%.4f'), 'alpha', num2str(alpha), '.eps'];

% Save the figure as an EPS file
print(gcf, filename, '-depsc');

% Close the figure window
% close(gcf);