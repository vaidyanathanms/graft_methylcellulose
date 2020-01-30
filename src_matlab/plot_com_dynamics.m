%% For plotting COM-time

clc;
clear;
close all;
format long;

%% Color Data

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',orange,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input Data

num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8'};
sigma_arr     = {'0.01','0.1','0.2','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];

fprintf('Plotting com distance as a function of time \n');

for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
    bb_mw  = mw_bb_arr(bb_cnt);
    
    for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
        gr_mw = mw_graft_arr(gr_cnt);
        
        for eps_cnt = 1:length(eps_arr) %epsilon loop
            eps_val = str2double(eps_arr{eps_cnt});
            
            hcom = figure;
            hold on
            box on
            set(gca,'FontSize',20)
            xlabel('Snapshot Number','FontSize',20,'Interpreter','Latex')
            ylabel('Mean Distance betweeen COM','FontSize',20,'Interpreter','Latex')
            
            for sig_cnt = 1:length(sigma_arr) %sigma loop
                
                sig_val = str2double(sigma_arr{sig_cnt});
                dirname = sprintf('../../sim_results/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
                    bb_mw,gr_mw,polydens,num_bb_chains);
                
                if ~exist(dirname,'dir')
                    fprintf('%s does not exist\n',dirname);
                    close all;
                    continue
                end
                    
                fprintf('Analyzing for bbMW/grMW/eps/sig: %g\t%g\t%g\t%g\n',...
                    bb_mw,gr_mw,sig_val,eps_val);
                
                com_prefix = sprintf('composallappend_%s_%s_%d.dat',...
                    eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                
                com_fylename = strcat(dirname,'/',com_prefix);
                
                if exist(com_fylename,'file') ~= 2
                    fprintf('%s does not exist\n',com_fylename);
                    continue;
                elseif struct(dir(com_fylename)).bytes == 0
                    fprintf('Empty file: %s \n',com_fylename);
                    continue;
                else
                    com_allvals = importdata(com_fylename);
                    lfyle = length(com_allvals.data(:,1));
                    plot(com_allvals.data(:,5),'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                    legendinfo{sig_cnt} = ['$\Sigma$: ' num2str(sig_val)];
                end
                
            end
            
            legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','NorthEast')
            saveas(hcom,sprintf('./../../all_figures/fig_dcomvstime_bbMW_%d_gMW_%d_nch_%d_eps_%s.png',...
                bb_mw,gr_mw,num_bb_chains,eps_arr{eps_cnt}))
            
        end
        
    end
    
end
