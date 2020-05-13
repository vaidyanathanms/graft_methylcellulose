%% To extract timestep from the input/job files

clc;
clear;
close all;
format long;

%% Input Data

num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8';'1.0';'1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [1];
from_backup   = 0; % If it needs to be read from backup drive - CHECK PATH


src_folder   = pwd;

for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    
    if from_backup == 1
        fprintf('Analyzing data from backup drive..\n');
    else
        fprintf('Analyzing data from "restart_files" folder ..\n');
    end
    
    for bb_cnt = 1:length(mw_bb_arr) % backbone MW loop  
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            for sig_cnt = 1:length(sigma_arr) %sigma loop
                sig_val = str2double(sigma_arr{sig_cnt});
                
                for eps_cnt = 1:length(eps_arr) %epsilon loop
                    eps_val = str2double(eps_arr{eps_cnt});
                    
                    if from_backup ~= 1
                        
                        simdirname = sprintf('../../restart_files/config_%d/restart_bbMW_%d_ngMW_%d_rho_%s_nch_%d/graftperc_%g/epsval_%s',...
                            config,bb_mw,gr_mw,polydens,num_bb_chains,sig_val,eps_arr{eps_cnt});
                    
                        if ~exist(simdirname,'dir') % Find if directory exists
                            fprintf('%s does not exist\n',simdirname);
                            continue
                        end
                        
                    else
                        
                        simdirname = sprintf(['H:/Trajectory_Backups/CG_GraftMethylCellulose/grafts_MC/newconfig_expandedlist/' ...
                            'config_%d/nchains_%d/backbonewtperc_%s/n_graft_%d/graftperc_%s/epsval_%s'],...
                            config,num_bb_chains,polydens,gr_mw,sigma_arr{sig_cnt},eps_arr{eps_cnt});

                        if ~exist(simdirname,'dir') % Find if directory exists
                            fprintf('%s does not exist\n',simdirname);
                            continue
                        end
                        
                    end    
                        
                    fprintf('Analyzing config/grMW/eps/sig: %d\t%g\t%g\t%g\n',config,gr_mw,eps_val,sig_val);
                 
                    job_prefix = sprintf('job*');
                    job_fylelist = dir(strcat(simdirname,'/',job_prefix));
                    
                    if min(size(job_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s in %s\n',job_prefix,simdirname);
                        break;
                    end
                    outfname = sprintf('../../timedata/out_timefiles/detailedtimevals_conf_%d_nch_%d_grmw_%d_eps_%s_sig_%s.dat',...
                        config,num_bb_chains,gr_mw,eps_arr{eps_cnt},sigma_arr{sig_cnt});
                    outf2name = sprintf('../../timedata/out_timefiles/summarytimevals_conf_%d_nch_%d_grmw_%d_eps_%s_sig_%s.dat',...
                        config,num_bb_chains,gr_mw,eps_arr{eps_cnt},sigma_arr{sig_cnt});
                    write_timevals(simdirname,job_fylelist,outfname,outf2name)
                    cd(src_folder);
                    
                end %end of epsilon loop
                
            end %end of sigma loop
            
        end % end of graft MW loop
        
    end %end of bb MW loop
    
end %end of config loop
