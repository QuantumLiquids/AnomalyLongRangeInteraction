% L_list = 6:2:28;
L_list = [6:19,20,22,24,26];
alpha = 1;
theta = 0.6;
marker_color1 = [19, 103, 131] / 256;
marker_color2 = [255, 158, 2] / 256;

gaps = zeros(1, numel(L_list));
mom_values = cell(1, numel(L_list));
for i = 1: numel(L_list)
    L = L_list(i);
    filename = ['EnergyN', num2str(L), 'theta', num2str(theta), 'alpha', num2str(alpha), '.txt'];
    % filename = ['EnergyN', num2str(L), 'alpha', num2str(alpha), '.txt'];
    % filename = ['EnergyN', num2str(L), '.txt'];
    % filename = ['EnergyN', num2str(L), 'theta', num2str(theta) '.txt'];
    energy_data = importdata(filename);
    energy_data_half = energy_data(1:L/2+1,:);
    E_list = energy_data_half(:);
    e012 = mink(E_list, 3);
    gaps(i) = e012(3) - e012(1);

    [row, col]=find(energy_data == e012(3));

    % E_list = data;
    % E0 = min(min(E_list));
    % gaps(i) = min(E_list(2,:))-E0;
    
    mom_values{i} =  ['$',char(sym((row - 1)/(L/2))),'\pi$'];
    fprintf('L : %i, ky_int : %i\n', L, row - 1 );
end

h = plot(1./L_list, gaps, 'o', 'MarkerSize', 10);hold on;
for i = 1:numel(L_list)
    L = L_list(i);
    T = text(1/L, gaps(i), mom_values{i}, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 18,'Interpreter','latex');
end


set(gca, 'fontsize', 24);
set(gca, 'linewidth', 1.5);
set(get(gca, 'Children'), 'linewidth', 2); % Set line width 1.5 pounds
xlabel('$1/L$', 'Interpreter', 'latex');
ylabel('gap $\Delta$', 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'FontSize', 24);
set(get(gca, 'YLabel'), 'FontSize', 24);
grid on;
xlim([0, inf]);
ylim([0, inf]);

% Perform linear fit
x = [1./L_list,0];
p = polyfit(1./L_list, gaps, 1);
fit_line = polyval(p, x);

hold on;
plot(x, fit_line, 'r-', 'LineWidth', 2);

% Calculate correlation coefficient (R-value)
R = corrcoef(1./L_list, gaps);
R_value = R(1, 2);

% Display R-value on the plot
text(0.6, 0.85, ['R = ', num2str(R_value)], 'Units', 'normalized', 'FontSize', 16);

legend('Data', 'Linear Fit', 'Location', 'best');