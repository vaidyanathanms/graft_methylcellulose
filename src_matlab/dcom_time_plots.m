%% To plot dcom for a given configuration vs time

clc;
clear;
close all;
format long;

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'k',orange,green,'m',brown,'b', gold};

%% Input Data

comflag       = [0,1,1,5];%[write,plot,XY-analysis_column] % Dont change this unless the files are changed
num_bb_chains = 2;
polydens      = '0.1';
eps_arr       = {'0.8'}%,'1.0','1.2'};
sigma_arr     = {'0.01','0.1','0.2','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [4];
maxstepcheck  = 4e7;
time_write    = 1; %to write time data from job files
dcom_0        = 2; %value for zero sigma (averaged value)

%% Main analysis - Write consolidated data

for conf_cnt = 1:length(config_arr) %config loop (replicate trials)
    
    config = config_arr(conf_cnt);
    fprintf('Configuration under analysis: %d\n', config);
    
    for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            fw_com = fopen(sprintf('../../dcom_data/maxsimtime/maxtime_config_%d_gMW_%d_nch_%d.dat',...
                config,gr_mw,num_bb_chains),'w');
            fprintf(fw_com,'%s\t %s\t %s\t %s\n','config','epsilon','sigma','maxtime');
            
            % Check main directory exists or not
            dirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
                config,bb_mw,gr_mw,polydens,num_bb_chains);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            for eps_cnt = 1:length(eps_arr) %sigma loop
                eps_val = str2double(eps_arr{eps_cnt});
                
                % Plot Data
                h_tcom = figure;
                hold on
                box on
                set( gcf, 'Units', 'normalized', 'Position', [0.1,0.1,0.8,0.8] ) ;
                hs = zeros(4);
                axis tight
                axis square
                x_min_max_times = zeros(4,2); %for setting x limits of plot
                y_min_max_times = zeros(4,2); %for setting y limits of plot
                
                for sig_cnt = 1:length(sigma_arr) %epsilon loop
                    sig_val = str2double(sigma_arr{sig_cnt});
                    
                    if time_write == 1 % begin time write
                        
                        com_prefix = sprintf('composallappend_%s_%s_%d.dat',...
                            eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                        com_fylename = strcat(dirname,'/',com_prefix);
                        
                        if exist(com_fylename,'file') ~= 2
                            fprintf('%s does not exist\n',com_fylename);
                            continue;
                        elseif struct(dir(com_fylename)).bytes == 0
                            fprintf('Empty file: %s \n',com_fylename);
                            continue;
                        end
                        
                        %load dcom data and sort data so that time jumbling is taken care of
                        xcol = comflag(1,3); ycol = comflag(1,4);
                        com_allvals = importdata(com_fylename);
                        com_allvals_unsrt = [com_allvals.data(:,xcol) com_allvals.data(:,ycol)];
                        com_allvals_srted = sortrows(com_allvals_unsrt,xcol);
                        clear com_allvals_unsrt
                        
                        % Compare with time-files and write dcom-time output to dcom-time directory
                        time_fyle = sprintf('../../timedata/out_timefiles/summarytimevals_conf_%d_nch_%d_grmw_%d_eps_%s_sig_%s.dat',...
                            config,num_bb_chains,gr_mw,eps_arr{eps_cnt},sigma_arr{sig_cnt});
                        
                        if exist(time_fyle,'file') ~= 2
                            fprintf('%s does not exist\n',time_fyle);
                            continue;
                        elseif struct(dir(time_fyle)).bytes == 0
                            fprintf('Empty file: %s \n',time_fyle);
                            continue;
                        end
                        
                        ftime_in  = importdata(time_fyle);
                        time_data = ftime_in.data;
                        
                        write_dcom_timevals(com_allvals_srted(:,1),com_allvals_srted(:,2),time_data,config,gr_mw,sig_val,eps_val)
                        
                    end %end time write
                    
                    dcom_fyle = importdata(sprintf('../../dcom_data/outdata/dcomtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                        config,gr_mw,sig_val,eps_val));
                    
                    if max(dcom_fyle.data(:,1)) < maxstepcheck
                        fprintf('WARNING: Not equilibrated \n');
                        fprintf('Config: %d/ eps: %g/ sig: %g/ maxtime: %g \n',config,eps_val,sig_val,max(com_allvals_srted(:,1)));
                        fprintf(fw_com,'%d\t%g\t%g\t%g\t%s\n',config,eps_val,sig_val,max(com_allvals_srted(:,1)),'WARNING: not equilibrated');
                    else
                        fprintf(fw_com,'%d\t%g\t%g\t%g\n',config,eps_val,sig_val,max(com_allvals_srted(:,1)));
                    end
                    
                    hs(sig_cnt) = subplot(4,1,sig_cnt);
                    axis square
                    plot(dcom_fyle.data(:,2),dcom_fyle.data(:,3),'Color',pclr{sig_cnt},...
                        'LineStyle','-','LineWidth',3)
                    
                    x_min_max_times(sig_cnt,1) = min(dcom_fyle.data(:,2));
                    x_min_max_times(sig_cnt,2) = max(dcom_fyle.data(:,2));
                    
                    y_min_max_times(sig_cnt,1) = min(dcom_fyle.data(:,3));
                    y_min_max_times(sig_cnt,2) = max(dcom_fyle.data(:,3));
                    
                    
                    lgd = legend(['\Sigma: ' num2str(sig_val)]);
                    lgd.FontSize = 14;
                    legend boxoff
                    
                end % end for each sigma value
                
                
                
                p1 = get(hs(1),'Position');
                p2 = get(hs(2),'Position');
                p3 = get(hs(3),'Position');
                p4 = get(hs(4),'Position');
                

                
                set(hs(1),'XTick',[],'YTick',[10 20 30]);
                set(hs(2),'XTick',[],'YTick',[10 20 30]);               
                set(hs(3),'XTick',[],'YTick',[10 20 30]);
                set(hs(4),'YTick',[10 20 30]);

                xlim(hs(1),[x_min_max_times(4,1) x_min_max_times(4,2)]);
                xlim(hs(2),[x_min_max_times(4,1) x_min_max_times(4,2)]);
                xlim(hs(3),[x_min_max_times(4,1) x_min_max_times(4,2)]);

                ymax = y_min_max_times(4,2)+0.1*y_min_max_times(4,2);
                ylim(hs(1),[y_min_max_times(4,1) ymax]);
                ylim(hs(2),[y_min_max_times(4,1) ymax]);
                ylim(hs(3),[y_min_max_times(4,1) ymax]);

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
                
                set(hs(1),'pos', p1,'FontSize',18);
                set(hs(2),'pos', p2,'FontSize',18);
                set(hs(3),'pos', p3,'FontSize',18);
                set(hs(4),'pos', p4,'FontSize',18);
                xlabel(hs(4),'Time ($\tau$)','FontSize',20,'Interpreter','Latex')
                ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
                
                saveas(h_tcom,sprintf('./../../all_figures/fig_dcomtime_conf_%d_gMW_%d_eps_%g.png',...
                    config,gr_mw,eps_val))
                
            end % end eps loop
            
        end % end graft loop
        
    end % end bb_mw loop
    
    fclose(fw_com);
    
end % end config loop

