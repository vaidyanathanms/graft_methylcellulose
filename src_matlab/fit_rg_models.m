clc;
clear;
close all;
format long;

%% Color Data

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Inputs

num_bb_chains = 1;
polydens      = '0.5';
eps_arr      = [0.8]%;1.0;1.2];
sigma_arr    = [0.01;0.05;0.1;0.15;0.2;0.3];
mw_graft_arr = [25]; %array should be of size 1
mw_bb_arr    = [800]; %array should be of size 1

%% Main analysis

% Load data into arrays
for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
    bb_mw  = mw_bb_arr(bb_cnt);
    
    for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
        gr_mw = mw_graft_arr(gr_cnt);
        
        rg_fylname = sprintf('../../outfiles/rgavg_bbMW_%d_gMW_%d_nch_%d.dat',...
            bb_mw,gr_mw,num_bb_chains);
        
        if exist(rg_fylname,'file') ~= 2
            fprintf('%s does not exist/empty file\n',rg_fylename);
            continue;
        end
        
        rg_all_data = importdata(rg_fylname);
        rgavg_sim   = extract_vals_sig_eps(rg_all_data,1,2,3,eps_arr,sigma_arr);
        rgsqavg_sim = rgavg_sim.^2;
        
        clear rg_all_data rg_fylname
        
        pers_fylname = sprintf('../../outfiles/pers_bbMW_%d_gMW_%d_nch_%d.dat',...
            bb_mw,gr_mw,num_bb_chains);
        
        if exist(pers_fylname,'file') ~= 2
            fprintf('%s does not exist/empty file\n',rg_fylename);
            continue;
        end
        
        pers_all_data = importdata(pers_fylname);
        pers_vals = extract_vals_sig_eps(pers_all_data,1,2,4,eps_arr,sigma_arr);
        clear pers_all_data pers_fylname
        
        term1 = 1/3*bb_mw*pers_vals;
        term2 = pers_vals.^2;
        term3 = 2*(pers_vals.^4./bb_mw^2).*(exp(-bb_mw./pers_vals)-1);
        term4 = 2*pers_vals.^3/bb_mw;
        rgsq_theory = term1 - term2 + term3 + term4;
        
        %Plot all data
        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
        ylabel('$R_{g}^2$','FontSize',20,'Interpreter','Latex')
        for ecnt = 1:length(eps_arr)
            plot(sigma_arr, rgsqavg_sim(ecnt,:),'MarkerSize', 10,...
                'LineStyle','none','MarkerFaceColor', pclr{ecnt},...
                'Marker', msty{ecnt},'MarkerEdgeColor',pclr{ecnt});
            plot(sigma_arr, rgsq_theory(ecnt,:),'MarkerSize', 10, ...
                'LineStyle','none','MarkerFaceColor', 'none',...
                'Marker', msty{ecnt},'MarkerEdgeColor',pclr{ecnt})
%             plot(sigma_arr,term1(ecnt,:),'r')
%             plot(sigma_arr,term2(ecnt,:),'b')
%             plot(sigma_arr,term3(ecnt,:),'g')
%             plot(sigma_arr,term4(ecnt,:),'k')
            legendinfo{2*ecnt-1} = ['Sim: $\epsilon_{pg}$ = ' num2str(eps_arr(ecnt))];
            legendinfo{2*ecnt} = ['Theory: $\epsilon_{pg}$ = ' num2str(eps_arr(ecnt))];
        end

        legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','Best')
        legend boxoff
        ylim([3 1.1*max(max(rgsq_theory))]);
        set(gca,'yscale','log')
        saveas(h1,sprintf('./../../all_figures/rgcomp_%d_bbMW_%d',gr_mw,bb_mw),'png');

        
        h1 = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
        ylabel('$R_{g}^2$','FontSize',20,'Interpreter','Latex')
        for ecnt = 1:length(eps_arr)
            plot(sigma_arr, rgsqavg_sim(ecnt,:),'MarkerSize', 10,...
                'LineStyle','none','MarkerFaceColor', pclr{ecnt},...
                'Marker', msty{ecnt},'MarkerEdgeColor',pclr{ecnt});
            plot(sigma_arr, rgsq_theory(ecnt,:)/sqrt(mw_bb_arr(1)),'MarkerSize', 10, ...
                'LineStyle','none','MarkerFaceColor', 'none',...
                'Marker', msty{ecnt},'MarkerEdgeColor',pclr{ecnt})
%             plot(sigma_arr,term1(ecnt,:),'r')
%             plot(sigma_arr,term2(ecnt,:),'b')
%             plot(sigma_arr,term3(ecnt,:),'g')
%             plot(sigma_arr,term4(ecnt,:),'k')
            legendinfo{2*ecnt-1} = ['Sim: $\epsilon_{pg}$ = ' num2str(eps_arr(ecnt))];
            legendinfo{2*ecnt} = ['Theory (Scaled): $\epsilon_{pg}$ = ' num2str(eps_arr(ecnt))];
        end

        legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','Best')
        legend boxoff
        ylim([3 1.1*max(max(rgsq_theory))]);
        set(gca,'yscale','log')
        saveas(h1,sprintf('./../../all_figures/rgcompscalesqrtN_%d_bbMW_%d',gr_mw,bb_mw),'png');

    end
    
end






