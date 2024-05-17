L_list = [8,10, 12, 14, 16, 18, 20];
alpha = 1;

marker_color1 = [19, 103, 131] / 256;
marker_color2 = [255, 158, 2] / 256;

e0s = zeros(1, numel(L_list));
for i = 1: numel(L_list)
    filename = ['EnergyN', num2str(L_list(i)), 'alpha', num2str(alpha), '.txt'];
    data = importdata(filename);

    E_list = data(:);
    e0s(i) = mink(E_list, 1);
end

h = plot(1./L_list, e0s, 'o', 'MarkerSize', 10);
set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2); % Set line width 1.5 pounds
xlabel('$1/L$', 'Interpreter', 'latex');
ylabel('Ground state energy $E_0$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 24);
set(get(gca, 'YLabel'), 'FontSize', 24);
grid on;
xlim([0, inf]);
% ylim([0, inf]);

% Perform linear fit
x = [1./L_list,0];
p = polyfit(1./L_list, e0s, 1);
fit_line = polyval(p, x);

hold on;
plot(x, fit_line, 'r-', 'LineWidth', 2);

% Calculate correlation coefficient (R-value)
R = corrcoef(1./L_list, e0s);
R_value = R(1, 2);

% Display R-value on the plot
text(0.6, 0.85, ['R = ', num2str(R_value)], 'Units', 'normalized', 'FontSize', 16);

legend('Data', 'Linear Fit', 'Location', 'best');