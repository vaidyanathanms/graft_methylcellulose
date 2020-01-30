% To append a set of files of the same format and write to a single file.
% Can be used to generate a single output file which can then be used to plot/analyze.
% Use analyze_static_properties.m to compute average and plot data obtained
% from this analysis.

clc;
close all;
clear;
format long;

%% Input Flags

rgflag    = [0,1,5,1]; %[flagid,columns_to_process-XY,num_headerlines];
shapeflag = [1,1]; %[flag,headerlines]
dcomflag  = [0,1]; %[flag,headerlines]
persflag  = [0,1,3,1]; %[flag,columns_to_process-XY,num_headerlines];

%% Input Data

num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8';'1.0';'1.2'};
sigma_arr     = {'0.0','0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
cutoff        = 0.002; %cutoff for persistence length
config_arr    = [1,2,4,5];

src_folder   = pwd;

for conf_cnt = 1:length(config_arr)
    
    config = config_arr(conf_cnt);
    
    for bb_cnt = 1:length(mw_bb_arr) % backbone MW loop
        
        bb_mw  = mw_bb_arr(bb_cnt);
        per_plot = zeros(length(sigma_arr),length(eps_arr),length(mw_graft_arr));
        
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
                    fprintf('Analyzing eps/sig/grMW/bbMW: %g\t%g\t%g\t%g\n',eps_val,sig_val,gr_mw,bb_mw);
                    
                    if rgflag(1,1) %% Find all files for rgflag
                        
                        rg_prefix = sprintf('rgavgall_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        rg_fylelist = dir(strcat(simdirname,'/',rg_prefix));
                        if size(rg_fylelist) == 0
                            fprintf('No files are found for %s\n',rg_prefix);
                            break;
                        end
                        
                        outfname = sprintf('rgavgallappend_%s_%s_%d.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        outfylename = strcat(simdirname,'/',outfname);
                        avgval = process_write_and_average(simdirname,rg_fylelist,rgflag,outfylename);
                        cd(src_folder);
                        
                    end %end rgflag
                    
                    distcomfile = -1;
                    if dcomflag(1,1) %%if distcom files are present
                        
                        distcomfile = 1;
                        com_prefix = sprintf('distcom_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        com_fylelist = dir(strcat(simdirname,'/',com_prefix));
                        if size(com_fylelist) == 0
                            fprintf('No files are found for %s\n',com_prefix);
                            break;
                        end
                        
                        outfname = sprintf('composallappend_%s_%s_%d.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        outfylename = strcat(simdirname,'/',outfname);
                        process_and_write_distCOM(simdirname,com_fylelist,dcomflag,outfylename,num_bb_chains);
                        cd(src_folder);
                        
                    end %end distCOM calculation
                    
                    if dcomflag(1,1) && distcomfile == -1 %COM calculation (compute iff distcom is not present)
                        
                        com_prefix = sprintf('compos_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        com_fylelist = dir(strcat(simdirname,'/',com_prefix));
                        if size(com_fylelist) == 0
                            fprintf('No files are found for %s\n',com_prefix);
                            break;
                        end
                        
                        outfname = sprintf('composallappend_%s_%s_%d.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        outfylename = strcat(simdirname,'/',outfname);
                        process_and_write_COM(simdirname,com_fylelist,dcomflag,outfylename,num_bb_chains);
                        cd(src_folder);
                        
                    end %end COM calculation
                    
                    if shapeflag(1,1) %Shape factor calculation
                        
                        shape_prefix = sprintf('eigMCavg_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        shape_fylelist = dir(strcat(simdirname,'/',shape_prefix));
                        
                        if size(shape_fylelist) == 0
                            fprintf('No files are found for %s\n',shape_prefix);
                            break;
                        end
                        outfname = sprintf('shapeallappend_%s_%s_%d.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        outfylename = strcat(simdirname,'/',outfname);
                        process_and_write_shapefactor(simdirname,shape_fylelist,shapeflag,outfylename,num_bb_chains);
                        cd(src_folder);
                        
                    end % end shape factor calculation
                    
                    if persflag(1,1) %begin persistence length calculation
                        
                        pers_prefix = sprintf('mainpersistautocf_%s_%s_%d_*.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        pers_fylelist = dir(strcat(simdirname,'/',pers_prefix));
                        
                        if size(pers_fylelist) == 0
                            fprintf('No files are found for %s\n',pers_prefix);
                            break;
                        end
                        
                        outfname = sprintf('mainpersistautocf_%s_%s_%d.dat',eps_arr{eps_cnt},...
                            sigma_arr{sig_cnt},bb_mw);
                        outfylename = strcat(simdirname,'/',outfname);
                        process_and_write_nth_file(simdirname,pers_fylelist,persflag,outfylename);
                        cd(src_folder);
                        
                    end %end persistence length calculation
                    
                end %end of epsilon loop
                
            end %end of sigma loop
            
        end % end of graft MW loop
        
    end %end of bb MW loop
  
end %end of config loop
