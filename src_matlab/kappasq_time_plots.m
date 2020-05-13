%% Plot kappasq-time, kappasq(t) autocorrelation, <kappasq(t)-mean(kappasq(t))> autocorrelation
clc;
clear;
close all;
format long

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {orange,green,'m','k',brown,'b','r',gold};
lsty = {'-','--',':'};

%% Input data

kappatime_flags = 1; %kappasq-time
autocorr_flags  = 0; %kappasq(t) autocorrelation, <kappasq(t)-mean(kappasq(t))> autocorrelation

num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8'};
sigma_arr     = {'0.01','0.1','0.2','0.3'};%'0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [2];
src_folder    = pwd;

for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    
    for bb_cnt = 1:length(mw_bb_arr) % backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            for eps_cnt = 1:length(eps_arr) %epsilon loop
                eps_val = str2double(eps_arr{eps_cnt});
                
                if kappatime_flags % Begin kappa_time plots
                    fprintf('Plotting kappasq_time for config/gr_mw/epsval: %d\t%d\t%g\n',...
                        config,gr_mw,eps_val)
                    
                    hksqtime = figure;
                    hold on
                    box on
                    set( gcf, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ) ;
                    set(gca,'FontSize',16)
                    x_min_max_times = zeros(4,2); %for setting x limits of plot
                    y_min_max_times = zeros(4,2); %for setting y limits of plot
                    hs = zeros(4);
                    for sig_cnt = 1:length(sigma_arr) %sigma loop
                        sig_val = str2double(sigma_arr{sig_cnt});
                        
                        fylename = sprintf('../../autocorr/kappasq_time/kappasqtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                            config,gr_mw,sig_val,eps_val);
                        if exist(fylename,'file') ~= 2
                            fprintf('%s does not exist\n',time_fyle);
                            continue;
                        elseif struct(dir(fylename)).bytes == 0
                            fprintf('Empty file: %s \n',time_fyle);
                            continue;
                        end
                        kappasq = importdata(fylename);
                        
                        hs(sig_cnt) = subplot(4,1,sig_cnt);
                        axis square
                        
                        plot(kappasq.data(:,2),kappasq.data(:,3),'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                        
                        
                        x_min_max_times(sig_cnt,1) = min(kappasq.data(:,2));
                        x_min_max_times(sig_cnt,2) = max(kappasq.data(:,2));
                        
                        y_min_max_times(sig_cnt,1) = min(kappasq.data(:,3));
                        y_min_max_times(sig_cnt,2) = max(kappasq.data(:,3));
                        
                        lgd = legend(['\Sigma: ' num2str(sig_val)]);
                        lgd.FontSize = 14;
                        legend boxoff
                        clear kappasq
                        
                    end % End sigma loop
                    
                    p1 = get(hs(1),'Position');
                    p2 = get(hs(2),'Position');
                    p3 = get(hs(3),'Position');
                    p4 = get(hs(4),'Position');
                    set(hs(1),'XTickLabel',[],'YTick',[0.25 0.5 0.75]);
                    set(hs(2),'XTickLabel',[],'YTick',[0.25 0.5 0.75]);
                    set(hs(3),'XTickLabel',[],'YTick',[0.25 0.5 0.75]);
                    set(hs(4),'YTick',[0.25 0.5 0.75]);
                    
                    xlim(hs(1),[x_min_max_times(4,1) x_min_max_times(4,2)]);
                    xlim(hs(2),[x_min_max_times(4,1) x_min_max_times(4,2)]);
                    xlim(hs(3),[x_min_max_times(4,1) x_min_max_times(4,2)]);
                    
                    
                    p1(4) = 0.775/4;
                    p2(4) = 0.775/4;
                    p3(4) = 0.775/4;
                    p4(4) = 0.775/4;
                    
                    p1(3) = 0.5*0.775;
                    p2(3) = 0.5*0.775;
                    p3(3) = 0.5*0.775;
                    p4(3) = 0.5*0.775;
                    
                    p4(2) = 0.15;
                    p3(2) = p4(2) + p4(4);
                    p2(2) = p3(2) + p3(4);
                    p1(2) = p2(2) + p2(4);
                    
                    set(hs(1),'pos', p1,'FontSize',16);
                    set(hs(2),'pos', p2,'FontSize',16);
                    set(hs(3),'pos', p3,'FontSize',16);
                    set(hs(4),'pos', p4,'FontSize',16);
                    xlabel(hs(4),'Time ($\tau$)','FontSize',22,'Interpreter','Latex')
                    ylabel('$\kappa^2$','FontSize',22,'Interpreter','Latex')
                    
                    saveas(hksqtime,sprintf('./../../all_figures/fig_ksqtime_conf_%d_gMW_%d_eps_%g.png',...
                        config,gr_mw,eps_val))
                    
                end % End kappa_time plots for a given eps_val
                
                
                if autocorr_flags % Begin autocorrelation plots
                    
                    fprintf('Plotting autocorrelation for config/gr_mw/epsval: %d\t%d\t%g\n',...
                        config,gr_mw,eps_val)
                    
                    hkacorr = figure; % Raw autocorrelation
                    hold on
                    box on
                    set(gca,'FontSize',20)
                    xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
                    ylabel('$\langle R_{T}(\kappa^2(t)) \rangle$','FontSize',20,'Interpreter','Latex')
                    
                    for sig_cnt = 1:length(sigma_arr) %sigma loop
                        sig_val = str2double(sigma_arr{sig_cnt});
                        
                        fylename = sprintf('../../autocorr/outvals/autocorr_gMW_%d_config_%d_sig_%g_eps_%g.dat',...
                            gr_mw,config,sig_val,eps_val);
                        if exist(fylename,'file') ~= 2
                            fprintf('%s does not exist\n',time_fyle);
                            continue;
                        elseif struct(dir(fylename)).bytes == 0
                            fprintf('Empty file: %s \n',time_fyle);
                            continue;
                        end
                        kappasq = importdata(fylename);
                        %plot only until half the timepoints available
                        tcutval = 0.9*max(kappasq.data(:,2));
                        for tfind = 1:length(kappasq.data(:,1))
                            if kappasq.data(tfind,2) > tcutval
                                tcutindex = tfind;
                                break;
                            end
                        end
                        
                        plot(kappasq.data(1:tcutindex,2),kappasq.data(1:tcutindex,4),...
                            'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                        legendinfo{sig_cnt} = ['$\Sigma$: ' num2str(sig_val)];
                        clear kappasq
                        
                    end % End sigma loop
                    
                    legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
                    saveas(hkacorr,sprintf('./../../all_figures/fig_ksqtime_conf_%d_gMW_%d_eps_%g.png',...
                        config,gr_mw,eps_val))
                    clear legendinfo
                    
                end % End autocorrelation data
                
            end %end of eps loop
            
        end % end of graft MW loop
        
    end %end of bbMW loop
    
end % end of config

