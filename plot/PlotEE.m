L = 128;
P = 1;
num_theta = 30;
theta= 0.31415926535897932384626433832795;% 0.94247779607693797153879301498385;
% 0.26179938779914943653855361527329% 
J = -1.6;
Level = 0;
Db = 200;


if(theta ~= 0)
    omega0_str = num2str(round(cos(theta),4));
    omega1_str = num2str(round(sin(theta),4));
else
    omega0_str = '1.0';
    omega1_str = '0.0';
end

file_name = ['../data/eeLRIL',num2str(L), 'omega0', omega0_str, 'omega1', omega1_str, 'J', num2str(round(J,1), '%.1f'), 'Level', num2str(Level), 'P', num2str(P), 'D', num2str(Db)];
file_id = fopen(file_name,'r');
ee2 = fread(file_id, L-1, 'double');

plot(ee2, '-o'); hold on;
l_list = 1:L-1;

start_site = 4;

modelfun = @(b,x)(b(1)/6 * log(sin(pi*(2 .* x + 1)/2./(L+1))) + b(2) - b(3) * sin(pi/2*(2.*x+1))./sqrt(sin(pi*(2 * x + 1)/2./(L+1))) );
mdl = fitnlm(l_list(start_site:1:end-start_site+1),ee2(start_site:1:end-start_site+1),modelfun,[1,0.7,sqrt(pi/L)]);

c = mdl.Coefficients.Estimate(1);
En = mdl.Coefficients.Estimate(2);
% kF = mdl.Coefficients.Estimate(4);
b = mdl.Coefficients.Estimate;
sites = l_list(start_site):0.1:l_list(end-start_site+1);
plot(sites, modelfun(b,sites),'-');

% plot(log((L+1)/pi * sin(pi*(2*l_list+1)/2/(L+1))),ee2,'-o');hold on;
% fit_x=log((L+1)/pi * sin(pi*(2*l_list(start_site:2:end-start_site)+1)/2/(L+1)));
% fit_y=ee2(start_site:2:end-start_site);
% p = fit((fit_x'),(fit_y'),'poly1');
% plot(fit_x, p.p1*fit_x + p.p2,'-.');hold on;
% if c < 0
%     c = 0;
% end
% if c > 1
%     c= 1;
% end

% plot([0,pi],[0,0],'r',  'LineWidth',2);hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$x$','Interpreter','latex','FontName','Times New Roman');
ylabel('Entanglement Entropy','Interpreter','latex','FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

set(gca,'linewidth',1.5);
set(gcf,'position',[1000,1000,750,600]);
