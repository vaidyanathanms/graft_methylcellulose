%% Block averaging for errors

clc;
clear;
close all;
format long;

%% Color Data

green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17];brown = [0.2 0 0];
pclr = {'m',orange,green,'k',brown, gold,'b','r'};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Input data for block averaging
num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'1.0','1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [4]; % Only for plotting
rg_bare       = 9.78;   % mean Rg of bare systems - no graft

%% Flags
rgflag      = [0,0,1,2];%[write,plot,XY-analysis_column]
timewrite   = 1; % write Rg values to output file at regular intervals
blockflag   = 1; % block average

%% Sampling Input Data
nsamples  = 3000;   % Sample size to initially look at -- NOTE: increase if this is not enough
tminval   = 10000;  % min time after which the data is looked at
tsample   = 6;      % time difference between nearest samples
xdim_data = 2;      % dimension of time data
ydim_data = 3;      % dimension of rg data

tmaxval = tminval + nsamples*tsample;
tinvals = tminval:tsample:tmaxval;

%% Main analysis - Write consolidated data

for conf_cnt = 1:length(config_arr) %config loop (replicate trials)
    
    config = config_arr(conf_cnt);
    fprintf('Configuration under analysis: %d\n', config);
    
    for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            % Check main directory exists or not
            dirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
                config,bb_mw,gr_mw,polydens,num_bb_chains);
            
            if ~exist(dirname,'dir')
                fprintf('%s does not exist\n',dirname);
                continue
            end
            
            for eps_cnt = 1:length(eps_arr) %sigma loop
                eps_val = str2double(eps_arr{eps_cnt});
                
                if blockflag
                    errfig = figure;
                    hold on
                    box on
                    set(gca,'FontSize',16)
                    xlabel('Block Size','FontSize',16,'Interpreter','Latex')
                    ylabel('Error','FontSize',16,'Interpreter','Latex')
                end
                
                for sig_cnt = 1:length(sigma_arr) %epsilon loop
                    sig_val = str2double(sigma_arr{sig_cnt});
                    
                    fprintf('Analyzing for config/bbMW/grMW/eps/sig: %d\t%g\t%g\t%g\t%g\n',...
                        config,bb_mw,gr_mw,sig_val,eps_val);
                    
                    rg_prefix = sprintf('rgavgallappend_%s_%s_%d.dat',eps_arr{eps_cnt},...
                        sigma_arr{sig_cnt},bb_mw);
                    rg_fylename = strcat(dirname,'/',rg_prefix);
                    
                    
                    if exist(rg_fylename,'file') ~= 2
                        fprintf('%s does not exist/empty file\n',rg_fylename);
                        continue;
                    elseif struct(dir(rg_fylename)).bytes == 0
                        fprintf('Empty file: %s \n',rg_fylename);
                        continue;
                    else
                        %find avg rg and write to file
                        rg_allvals = importdata(rg_fylename);
                        xcol = rgflag(1,3); ycol = rgflag(1,4);
                        rgunsrtdata = [rg_allvals.data(:,xcol) rg_allvals.data(:,ycol)];
                        sortrg = sortrows(rgunsrtdata,1); %to account for the fact that the data maybe jumbled
                        lenrg  = length(sortrg(:,1));
                        
                        if timewrite
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
                            write_rg_timevals(sortrg(:,1), sortrg(:,2), time_data, config,gr_mw,sig_val,eps_val);
                            
                        end
                        
                        if blockflag
                            
                            fprintf('Sampling time series for config/gr_mw/epsval: %d\t%d\t%g\n',...
                                config,gr_mw,eps_val)
                            
                            rg_time_fyle = sprintf('../../rgtime_data/time_data/rgtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                                config,gr_mw,sig_val,eps_val);
                            if exist(rg_time_fyle,'file') ~= 2
                                fprintf('%s does not exist\n',time_fyle);
                                continue;
                            elseif struct(dir(rg_time_fyle)).bytes == 0
                                fprintf('Empty file: %s \n',time_fyle);
                                continue;
                            end
                            rgsq = importdata(rg_time_fyle);
                            
                            rgtimearr = sample_rg(rgsq.data, xdim_data, ydim_data, tinvals, tsample, config, gr_mw, sig_val, eps_val);
                            std_out   = compute_block_averages(rgtimearr(:,2));
                            plot(std_out(:,1),std_out(:,3),'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                            legendinfo{sig_cnt} = ['$\Sigma$: ' sigma_arr{sig_cnt}];
                            clear rgtimearr std_out
                            
                        end
                        
                    end
                    
                    clear rg_fylename rg_time_fyle
                    
                end % end sigma loop
                
                if blockflag
                    legend(legendinfo,'Interpreter','Latex','FontSize',12,'Location','Best');
                    saveas(errfig,sprintf('./../../all_figures/fig_errordata_conf_%d_gMW_%d_nch_%d_eps_%s.png',...
                        config,gr_mw,num_bb_chains,eps_arr{eps_cnt}));
                    clear legendinfo
                end
                
            end % end epsilon loop
            
        end %end graft MW loop
        
    end %end backbone MW loop
    
end %end config loop