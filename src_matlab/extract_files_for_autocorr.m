%% To extract files where trajectories are saved at a shorter interval
% Can be then used to compute autocorrelation

clc;
clear;
close all;
format long;

%% Input Data

dt_frames     = 100; %extract files with this step difference

num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8','1.0','1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [1,2,4,5];

src_folder    = pwd;
shapeflag     = [1,1]; %[flag,headerlines] DO NOT CHANGE THIS

%% Main analysis - Write consolidated data

for conf_cnt = 1:length(config_arr) %config loop (replicate trials)
    
    config = config_arr(conf_cnt);
    fprintf('Configuration under analysis: %d\n', config);
    
    for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            simdirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
                config,bb_mw,gr_mw,polydens,num_bb_chains);
            
            if ~exist(simdirname,'dir') % Find if directory exists
                fprintf('%s does not exist\n',simdirname);
                continue
            end
            
            for sig_cnt = 1:length(sigma_arr) %sigma loop
                sig_val = str2double(sigma_arr{sig_cnt});
                
                for eps_cnt = 1:length(eps_arr) %epsilon loop
                    eps_val = str2double(eps_arr{eps_cnt});
                    
                    fprintf('Analyzing short-time autocf eps/sig/grMW/bbMW: %g\t%g\t%g\t%g\n',eps_val,sig_val,gr_mw,bb_mw);
                    outfname = sprintf('../../autocorr/short_time_data/autocfshort_conf_%d_eps_%s_sig_%s.dat',...
                        config,eps_arr{eps_cnt},sigma_arr{sig_cnt});
                    
                    %Initiate output files and get file ID
                    shape_prefix = sprintf('eigMCavg_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                        sigma_arr{sig_cnt},bb_mw);
                    shape_fylelist = dir(strcat(simdirname,'/',shape_prefix));
                    
                    if min(size(shape_fylelist)) == 0
                        fprintf('No files/Empty files are found for %s\n',shape_prefix);
                        break;
                    end
                    
                    check_time_diff_and_write(simdirname,shape_fylelist,shapeflag,outfname,dt_frames);
                    cd(src_folder);
                    
                end
                
            end
            
        end
        
    end
    
end
                    
                    
                    