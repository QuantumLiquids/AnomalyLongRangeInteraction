L = 128;
P = 1;
num_theta = 30;
theta_list= [0:(pi/2)/num_theta:pi/2]; theta_list(end)=[];
J_list = [-2:0.1:1.9];
central_charge_list = zeros(numel(theta_list), numel(J_list));
Level = 0;
possible_Db = [2000,1000,600,200];

central_charge_filename = ['central_charge/centralchargeL', num2str(L), 'Level', num2str(Level), 'P', num2str(P),'.mat'];
if(0 && exist(central_charge_filename,'file'))
    load(central_charge_filename);
else
    for i = 1:numel(theta_list)
        theta = theta_list(i);
        if(theta ~= 0)
            omega0_str = num2str(round(cos(theta),4));
            omega1_str = num2str(round(sin(theta),4));
        else
            omega0_str = '1.0';
            omega1_str = '0.0';
        end
        for j = 1:numel(J_list)
            J = J_list(j);
            for k = 1:numel(possible_Db)
                Db = possible_Db(k);
                file_name = ['../data/eeLRIL',num2str(L), 'omega0', omega0_str, 'omega1', omega1_str, 'J', num2str(round(J,1), '%.1f'), 'Level', num2str(Level), 'P', num2str(P), 'D', num2str(Db)];
                if( exist(file_name,'file') )
                    break;
                end
            end
            file_id = fopen(file_name,'r');
            ee2 = fread(file_id,L-1, 'double');
            l_list = 1:L-1;

            start_site = 4;

            modelfun = @(b,x)(b(1)/6 * log(sin(pi*(2 .* x + 1)/2./(L+1))) + b(2) - b(3) * sin(pi/2*(2.*x+1))./sqrt(sin(pi*(2 * x + 1)/2./(L+1))) );
            mdl = fitnlm(l_list(start_site:1:end-start_site+1),ee2(start_site:1:end-start_site+1),modelfun,[1,0.7,sqrt(pi/L)]);

            c = mdl.Coefficients.Estimate(1);
            En = mdl.Coefficients.Estimate(2);
            % kF = mdl.Coefficients.Estimate(4);
            b = mdl.Coefficients.Estimate;
            sites = l_list(start_site):0.1:l_list(end-start_site+1);


            if c < -1 ||  c > 2
                fprintf('theta = %.4f, J = %.1f, c = %.2f, D = %i \n',theta, J, c, Db);
            end
            central_charge_list(i,j)=c;
        end
    end
    save(central_charge_filename, 'central_charge_list');
end
imagesc(J_list, theta_list, central_charge_list);hold on;
colorbar
set(gca,'YDir','normal');


% plot([0,pi],[0,0],'r',  'LineWidth',2);hold on;

set(gca,'fontsize',24);
set(gca,'linewidth',1.5);
%set(get(gca,'Children'),'linewidth',2); % Set line width 1.5 pounds
xlabel('$J$','Interpreter','latex','FontName','Times New Roman');
ylabel('$\theta$','Interpreter','latex','FontName','Times New Roman');
set(get(gca,'XLabel'),'FontSize',24);
set(get(gca,'YLabel'),'FontSize',24);

set(gca,'linewidth',1.5);
set(gcf,'position',[1000,1000,750,600]);
