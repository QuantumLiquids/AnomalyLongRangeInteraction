marker_color1 =  [019, 103, 131]/256;
marker_color2 =  [255,158,002] / 256;

data = importdata('data.txt');
E_list = data;
E0 = min(min(E_list));
gap1 = min(E_list(2,:))-E0;
E_list = E_list - E0;
E_list = E_list /gap1;

L = size(E_list,1);
E_list = [E_list(L/2+1:end,:); E_list(1:L/2,:)];

% Plot the E values for the selected k
h_list = plot(-L/2:L/2-1,E_list,'o', 'MarkerSize',10);
for h = h_list
    set(h, 'MarkerFaceColor', marker_color2);
    set(h, 'MarkerEdgeColor', marker_color2);
end
set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$S$','Interpreter','latex');
ylabel('$\Delta$','Interpreter','latex');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

grid on;
% title(['Energy Spectrum for k = ', num2str(k_index)]);