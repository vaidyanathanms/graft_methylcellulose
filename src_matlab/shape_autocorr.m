%% Compute autocorrelation of shape factors

clc;
clear;
close all;
format long

%% Input data

shapedata     = {'shapefactor',7}; % key for the shapefactor and the column
write_time    = 1; %1 - if time data has to be written; 0 - only autocorrelation
num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8'};%'0.8';'1.0';'1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [2];
src_folder   = pwd;
tcut_auto    = 1e3; %time after which mean is computed

for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    
    for bb_cnt = 1:length(mw_bb_arr) % backbone MW loop
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
            
            for sig_cnt = 1:length(sigma_arr) %sigma loop
                sig_val = str2double(sigma_arr{sig_cnt});
                
                for eps_cnt = 1:length(eps_arr) %epsilon loop
                    eps_val = str2double(eps_arr{eps_cnt});
                    
                    fprintf('Computing autocorrelation for config/gr_mw/sigval/epsval: %d\t%d\t%g\t%g\n',...
                        config,gr_mw,sig_val,eps_val)
                    
                    if write_time % begin writing time data for kappasq
                        %read kappa^2 data
                        shape_prefix = sprintf('shapeallappend_%s_%s_%d.dat',...
                            eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                        
                        kappasq_unsrt_fyle = strcat(dirname,'/',shape_prefix);
                        
                        if exist(kappasq_unsrt_fyle,'file') ~= 2
                            fprintf('%s does not exist\n',kappasq_unsrt_fyle);
                            continue;
                        elseif struct(dir(kappasq_unsrt_fyle)).bytes == 0
                            fprintf('Empty file: %s \n',kappasq_unsrt_fyle);
                            continue;
                        end
                        
                        shape_allvals = importdata(kappasq_unsrt_fyle);
                        % check whether the shapefactor column is same as given
                        % in input data
                        if strcmp(shape_allvals.colheaders{shapedata{2}},shapedata{1}) == 0
                            fprintf('Key for shape factor not found %s\n', shape_allvals.text)
                            continue
                        end
                        
                        % Sort data
                        shapeunsrt = shape_allvals.data;
                        shapesrt = sortrows(shapeunsrt,1);
                        
                        % Compare with time-files and write kappasq-time output
                        % to kappasq-time directory
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
                        
                        % inputs the kappasq values and the simulation time
                        % data. Outputs the kappasq data with duplicates
                        % removed. Will write all the data. Needs to be reread
                        % for autocorrelation calculation.
                        write_kappasq_timevals(shapesrt(:,1),shapesrt(:,shapedata{2}),time_data,config,gr_mw,sig_val,eps_val);
                        
                    end % End timewrite (kappasq data)
                    
                    
                    % Begin autocorrelation calculation
                    shape_fyle = sprintf('../../autocorr/kappasq_time/kappasqtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                        config,gr_mw,sig_val,eps_val);
                    if exist(shape_fyle,'file') ~= 2
                        fprintf('%s does not exist\n',time_fyle);
                        continue;
                    elseif struct(dir(shape_fyle)).bytes == 0
                        fprintf('Empty file: %s \n',time_fyle);
                        continue;
                    end
                    shape_fyle = importdata(shape_fyle);
                    kappasqout = shape_fyle.data(:,3);
                    
                    % Compute autocorrelation
                    shape_autocf = compute_autocorr(kappasqout,1);
                    
                    % write autocorrelation output
                    fw_com = fopen(sprintf('../../autocorr/outvals/autocorr_gMW_%d_config_%d_sig_%g_eps_%g.dat',...
                        gr_mw,config,sig_val,eps_val),'w');
                    fprintf(fw_com,'%s\t%s\t%s\t%s\n','Step','Time','Autocorr','Normautocorr');
                    fprintf(fw_com,'%g\t%g\t%g\t%g\n',[shape_fyle.data(:,1) shape_fyle.data(:,2) shape_autocf shape_autocf/shape_autocf(1)]');
                    fclose(fw_com);
                    
                end %end of epsilon loop
                
            end %end of sigma loop
            
        end % end of graft MW loop
        
    end %end of bb MW loop
    
end %end of config loop
