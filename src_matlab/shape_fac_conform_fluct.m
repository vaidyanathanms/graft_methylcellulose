%% To analyze shape factors and conformational fluctuations

clc;
clear;
close all;
format long;

%% Input Data

num_bb_chains = 2;
polydens      = '0.1';
eps_arr       = {'0.8'};%,'1.0','1.2'};
sigma_arr     = {'0.2'};%,'0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config        = 1;

cutoff       = 0.1; %cutoff for persistence length
tcut_frac    = 0.8; %cutoff time fraction to start calculating dcomavg
tcut_prop_avg = 10^7; %cutoff timestep beyond which avg values are computed


%% Main analysis - Write consolidated data

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
            
            for sig_cnt = 1:length(sigma_arr) %epsilon loop
                sig_val = str2double(sigma_arr{sig_cnt});
                
                fprintf('Analyzing for config/bbMW/grMW/eps/sig: %d\t%g\t%g\t%g\t%g\n',...
                    config,bb_mw,gr_mw,sig_val,eps_val);
                
                fprintf('Analyzing shape factor values \n');
                shape_prefix = sprintf('shapeallappend_%s_%s_%d.dat',...
                    eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                
                shape_fylename = strcat(dirname,'/',shape_prefix);
                
                if exist(shape_fylename,'file') ~= 2
                    fprintf('%s does not exist\n',shape_fylename);
                    continue;
                elseif struct(dir(shape_fylename)).bytes == 0
                    fprintf('Empty file: %s \n',shape_fylename);
                    continue;
                else
                    
                    shape_allvals = importdata(shape_fylename);
                    lfyle = length(shape_allvals.data(:,1));
                    shapeunsrt = shape_allvals.data;
                    shapesrt = sortrows(shapeunsrt,1); %sort for time jumbling if any
                    
                    %Create and plot histograms and kappa^2-time data
                    create_and_plot_fluctuation_histogram(shapesrt,num_bb_chains,lfyle,config,bb_mw,gr_mw,eps_val,sig_val)

                    %compute averages
                    lenkappa   = length(shapesrt(:,1));
                    tcut_index = find_general_cutoff_index(shapesrt(:,1),tcut_prop_avg);
                    avg_sphere = zeros(num_bb_chains,1);
                    avg_cylndr = zeros(num_bb_chains,1);
                    avg_shape  = zeros(num_bb_chains,1);
                    
                    ncntr = zeros(num_bb_chains,1);
                    
                    for cntr = tcut_index:lfyle
                        chid = shapesrt(cntr,2);
                        avg_sphere(chid,1) = avg_sphere(chid,1) + shapesrt(cntr,5);
                        avg_cylndr(chid,1) = avg_cylndr(chid,1) + shapesrt(cntr,6);
                        avg_shape(chid,1) = avg_shape(chid,1) + shapesrt(cntr,7);
                        ncntr(chid,1) = ncntr(chid,1) + 1;
                    end
                    ncntr = ncntr/num_bb_chains;
                    for chainID = 1:num_bb_chains
                        sph_mean = avg_sphere(chainID,1)/ncntr(chainID,1);
                        cyl_mean = avg_cylndr(chainID,1)/ncntr(chainID,1);
                        shp_mean = avg_shape(chainID,1)/ncntr(chainID,1);
                    end    
                end
                
                clear shape_fylename
                clear avg_sphere avg_cylndr avg_shape
                clear xcol ycol tcut_index
                
            end
        end
    end
end
                
                