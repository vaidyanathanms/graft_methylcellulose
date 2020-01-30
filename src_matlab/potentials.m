%% To compare potentials between different monomers

clc;
clear;
close all;
format long;


%% Color Data
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k',orange, gold,'r','b'};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs

eps_graft = [0.8;1.0;1.2]; sig_graft = 0.823916; rc_graft = 1.352094;

eps_mc = [1.278034;2.093957;1.932788;1.566994];
sig_mc = [1.198800;1.068700;1.300100;1.452100];
rc_mc  = [1.969500;1.747600;2.135900;2.385600];

rvals = 0.75:0.03:2.5;

%% Main calculation and plots

h1 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$r (\sigma)$','FontSize',20,'Interpreter','Latex')
ylabel('$U (k_B T)$','FontSize',20,'Interpreter','Latex');


for gcnt = 1:length(eps_graft)
    epsval = eps_graft(gcnt); sigval = sig_graft; 
    uatcutoff = 4*epsval*( ((sigval/rc_graft).^9) - ((sigval/rc_graft).^6) ); 
    uval = 4*epsval.*( ((sigval./rvals).^9) - ((sigval./rvals).^6) ) - uatcutoff; 
    for i = 1:length(rvals)
        if rvals(i) > rc_graft
            uval(i) = 0;
        end
    end
    plot(rvals,uval,'Color',pclr{gcnt},'LineWidth',2,'LineStyle',lsty{1})
    legendinfo{gcnt} = ['$\epsilon_{gg}$ = ' num2str(epsval)];
end


for mcnt = 1:length(eps_mc)
    epsval = eps_mc(mcnt); sigval = sig_mc(mcnt); rcval = rc_mc(mcnt);
    uatcutoff = 4*epsval*( ((sigval/rcval).^9) - ((sigval/rcval).^6) ); 
    uval = 4*epsval.*( ((sigval./rvals).^9) - ((sigval./rvals).^6) ) - uatcutoff; 
    for i = 1:length(rvals)
        if rvals(i) > rcval
            uval(i) = 0;
        end
    end
    plot(rvals,uval,'Color',pclr{mcnt+gcnt},'LineWidth',2,'LineStyle',lsty{2})
end

legendinfo{gcnt+1} = ['$C_0, \epsilon_{bb}$ = ' num2str(eps_mc(1))];
legendinfo{gcnt+2} = ['$C_1, \epsilon_{bb}$ = ' num2str(eps_mc(2))];
legendinfo{gcnt+3} = ['$C_2, \epsilon_{bb}$ = ' num2str(eps_mc(3))];
legendinfo{gcnt+4} = ['$C_3, \epsilon_{bb}$ = ' num2str(eps_mc(4))];

zerovals = zeros(length(rvals),1);
plot(rvals,zerovals,'--b','LineWidth',1)

legend(legendinfo,'Interpreter','Latex','FontSize',14,'Location','Best')
legend boxoff
ylim([-1.5 2]);
xlim([min(rvals) max(rvals)]);

