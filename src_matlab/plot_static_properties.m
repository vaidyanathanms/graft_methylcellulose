%% To analyze static properties

clc;
clear;
close all;
format long;

%% Color Data

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',brown,green,'k','m', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input Flags

rgflag   = 1;
persflag = 1;
comp_MW  = 0;
%% Input Data

num_bb_chains = 1;
polydens      = '0.5';
eps_arr      = [0.8;1.0;1.2];
sigma_arr    = [0.01;0.05;0.10;0.15;0.2;0.25;0.30];
mw_graft_arr = [50];
mw_bb_arr    = [1000];

cutoff       = 0.001; %cutoff for persistence length

if persflag
    
    for bb_cnt = 1:length(mw_bb_arr)
        bb_mw  = mw_bb_arr(bb_cnt);
        per_plot = zeros(length(sigma_arr),length(eps_arr),length(mw_graft_arr));
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            
            gr_mw = mw_graft_arr(gr_cnt);            

            h1 = figure;
            hold on
            box on
            set(gca,'FontSize',16)
            xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
            ylabel('$l_{p}$','FontSize',20,'Interpreter','Latex')
                    
            per_fyle = sprintf('../../outfiles/pers_bbMW_%d_gMW_%d_nch_%d.dat',...
                bb_mw,gr_mw,num_bb_chains);
            
            if exist(per_fyle,'file') ~= 2 
                fprintf('%s does not exist\n', per_fyle);
                continue;
            end
            
            fprintf('Analyzing %s\n',per_fyle);
            
            alldata  = importdata(per_fyle);
            sig_data = alldata.data(:,1); %Column 1 is sigma
            eps_data = alldata.data(:,2); %Column 2 is epsilon
            len_data = length(sig_data(:,1));
            
            for datacnt = 1:len_data
                
                f_flag = -1;
                for f_dat = 1:length(sigma_arr)
                
                    if sig_data(datacnt,1) == sigma_arr(f_dat,1)
                        
                        sig_ind = f_dat; f_flag = 1;
                        break;
                        
                    end
                    
                end
                
                if f_flag == -1 
                    
                    fprintf('Did not find corresponding sigma value for %g',...
                        sig_data(datacnt,1));
                    continue;
                    
                end
                
                f_flag = -1;
                
                for f_dat = 1:length(eps_arr)
                    
                    if eps_data(datacnt,1) == eps_arr(f_dat,1)
                        
                        eps_ind = f_dat; f_flag = 1;
                        break;

                    end
                    
                end
                
                if f_flag == -1 
                    
                    fprintf('Did not find corresponding sigma value for %g',...
                        sig_data(datacnt,1));
                    continue;
                    
                end
                
                per_plot(sig_ind,eps_ind,gr_cnt) = alldata.data(datacnt,4); %4th column is persistence length
                
            end
           
            % Plot as a function of epsilon for each MW
            for epscnt = 1:length(eps_arr)
                
                plot(sigma_arr, per_plot(:,epscnt,gr_cnt),'MarkerSize', 10, ...
                    'LineStyle','none','MarkerFaceColor', pclr{epscnt},...
                    'Marker', msty{epscnt},'MarkerEdgeColor', pclr{epscnt})
            
                legendinfo{epscnt} = ['$\epsilon_{pg}$ = ' num2str(eps_arr(epscnt))];
    
            end
            
            legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','Best')
            legend boxoff
            ylim([3 1.1*max(max(max(per_plot)))]);
            saveas(h1,sprintf('./../../all_figures/gMW_%d_bbMW_%d',gr_mw,bb_mw),'png');
            
        end
        clear legendinfo
        if comp_MW
            %Compare between gr_MW
            for epscnt = 1:length(eps_arr)
                
                h1 = figure;
                hold on
                box on
                set(gca,'FontSize',16)
                xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
                ylabel('$l_{p}$','FontSize',20,'Interpreter','Latex')
                
                for mwcnt = 1:length(mw_graft_arr)
                    
                    plot(sigma_arr, per_plot(:,epscnt,mwcnt),'MarkerSize', 10, ...
                        'LineStyle','none','MarkerFaceColor', pclr{mwcnt},...
                        'Marker', msty{mwcnt},'MarkerEdgeColor', pclr{mwcnt})
                    
                    legendinfo{mwcnt} = ['MW$_{graft}$ = ' num2str(mw_graft_arr(mwcnt))];
                    
                end
                
                legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','NorthWest')
                legend boxoff
                ylim([3 1.1*max(max(max(per_plot)))]);
                saveas(h1,sprintf('./../../all_figures/eps_%g_bbMW_%d',eps_arr(epscnt),bb_mw),'png');
                
                clear legendinfo
                
            end
            
        end
            
    end
    
end


if rgflag
    
     for bb_cnt = 1:length(mw_bb_arr)
        bb_mw  = mw_bb_arr(bb_cnt);
        per_plot = zeros(length(sigma_arr),length(eps_arr),length(mw_graft_arr));
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            
            gr_mw = mw_graft_arr(gr_cnt);            

            h1 = figure;
            hold on
            box on
            set(gca,'FontSize',16)
            xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
            ylabel('$\langle R_{g}^2 \rangle^{1/2}$','FontSize',20,'Interpreter','Latex')
                    
            per_fyle = sprintf('../../outfiles/rgavg_bbMW_%d_gMW_%d_nch_%d.dat',...
                bb_mw,gr_mw,num_bb_chains);
            
            if exist(per_fyle,'file') ~= 2 
                fprintf('%s does not exist\n', per_fyle);
                continue;
            end
            
            fprintf('Analyzing %s\n',per_fyle);
            
            alldata  = importdata(per_fyle);
            sig_data = alldata.data(:,1); %Column 1 is sigma
            eps_data = alldata.data(:,2); %Column 2 is epsilon
            len_data = length(sig_data(:,1));
            
            for datacnt = 1:len_data
                
                f_flag = -1;
                for f_dat = 1:length(sigma_arr)
                
                    if sig_data(datacnt,1) == sigma_arr(f_dat,1)
                        
                        sig_ind = f_dat; f_flag = 1;
                        break;
                        
                    end
                    
                end
                
                if f_flag == -1 
                    
                    fprintf('Did not find corresponding sigma value for %g',...
                        sig_data(datacnt,1));
                    continue;
                    
                end
                
                f_flag = -1;
                
                for f_dat = 1:length(eps_arr)
                    
                    if eps_data(datacnt,1) == eps_arr(f_dat,1)
                        
                        eps_ind = f_dat; f_flag = 1;
                        break;

                    end
                    
                end
                
                if f_flag == -1 
                    
                    fprintf('Did not find corresponding sigma value for %g',...
                        sig_data(datacnt,1));
                    continue;
                    
                end
                
                per_plot(sig_ind,eps_ind,gr_cnt) = alldata.data(datacnt,3); %3rd column is Rg
                
            end
           
            % Plot as a function of epsilon for each MW
            for epscnt = 1:length(eps_arr)
                
                plot(sigma_arr, per_plot(:,epscnt,gr_cnt),'MarkerSize', 10, ...
                    'LineStyle','none','MarkerFaceColor', pclr{epscnt},...
                    'Marker', msty{epscnt},'MarkerEdgeColor', pclr{epscnt})
            
                legendinfo{epscnt} = ['$\epsilon_{pg}$ = ' num2str(eps_arr(epscnt))];
    
            end
            
            legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','Best')
            legend boxoff
            ylim([3 1.1*max(max(max(per_plot)))]);
            saveas(h1,sprintf('./../../all_figures/gMW_%d_bbMW_%d',gr_mw,bb_mw),'png');
            
        end
        clear legendinfo
        if comp_MW
            %Compare between gr_MW
            for epscnt = 1:length(eps_arr)
                
                h1 = figure;
                hold on
                box on
                set(gca,'FontSize',16)
                xlabel('$\sigma$','FontSize',20,'Interpreter','Latex')
                ylabel('$\langle R_{g}^2 \rangle^{1/2}$','FontSize',20,'Interpreter','Latex')
                
                for mwcnt = 1:length(mw_graft_arr)
                    
                    plot(sigma_arr, per_plot(:,epscnt,mwcnt),'MarkerSize', 10, ...
                        'LineStyle','none','MarkerFaceColor', pclr{mwcnt},...
                        'Marker', msty{mwcnt},'MarkerEdgeColor', pclr{mwcnt})
                    
                    legendinfo{mwcnt} = ['MW$_{graft}$ = ' num2str(mw_graft_arr(mwcnt))];
                    
                end
                
                legend(legendinfo,'Interpreter','Latex','FontSize',24,'Location','NorthWest')
                legend boxoff
                ylim([3 1.1*max(max(max(per_plot)))]);
                saveas(h1,sprintf('./../../all_figures/eps_%g_bbMW_%d',eps_arr(epscnt),bb_mw),'png');
                
                clear legendinfo
                
            end
            
        end
        
    end
    
end
    
    
