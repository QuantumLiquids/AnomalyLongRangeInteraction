L_list = [8:2:26];
alpha = 1;
theta = 0.6;
marker_color1 =  [019, 103, 131]/256;
marker_color2 =  [255,158,002] / 256;
marker_color3 =  [251,056,071] / 256;
marker_color4 =  [131,064,028] / 256;

gaps = zeros(1, numel(L_list));
mom_values = cell(1, numel(L_list));
for i = 1: numel(L_list)
    L = L_list(i);
    filename = ['EnergyN', num2str(L), 'theta', num2str(theta), 'alpha', num2str(alpha), '.txt'];
    % filename = ['EnergyN', num2str(L), 'alpha', num2str(alpha), '.txt'];
    % filename = ['EnergyN', num2str(L), '.txt']; %local hamiltonian, theta=pi/4
    % filename = ['EnergyN', num2str(L), 'theta', num2str(theta) '.txt'];
    energy_data = importdata(filename);
    if(mod(L,2)==0)
        energy_data_half = energy_data(1:L/2+1,:);
    else
        energy_data_half = energy_data(1:(L+1)/2,:);
    end
    E_list = energy_data_half(:);
    e01 = mink(E_list, 2);
    gaps(i) = e01(2) - e01(1);
    % gaps(i)= min(energy_data_half(2,:))-min(energy_data_half(1,:));

    [row, col]=find(energy_data == e01(2));

    % E_list = data;
    % E0 = min(min(E_list));
    % gaps(i) = min(E_list(2,:))-E0;
    if(mod(L,2)==0)
        mom_values{i} =  ['$',char(sym((row - 1)/(L/2))),'\pi$'];
    else
        mom_values{i} =  ['$',char(sym(2*(row - 1)/L)),'\pi$'];
    end
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