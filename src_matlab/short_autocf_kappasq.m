%% Compute autocorrelation of shape factors at short times

clc;
clear;
close all;
format long

%% Input data

shapedata     = {'shapefactor',7}; % key for the shapefactor and the column
dtvals        = 0.0003; %for now, all dt are same. If not it needs to be given as an array of size eps_arr*sigma_arr
num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8';'1.0';'1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [5];
src_folder   = pwd;

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
                    
                    % Begin autocorrelation calculation
                    shape_fyle = sprintf('../../autocorr/short_time_data/autocfshort_conf_%d_eps_%s_sig_%s.dat',...
                        config,eps_arr{eps_cnt},sigma_arr{sig_cnt});
                    if exist(shape_fyle,'file') ~= 2
                        fprintf('%s does not exist\n',shape_fyle);
                        continue;
                    elseif struct(dir(shape_fyle)).bytes == 0
                        fprintf('Empty file: %s \n',shape_fyle);
                        continue;
                    end
                    
                    shape_fyle = importdata(shape_fyle);
                    kappasqout = shape_fyle.data(:,7);
                    time_arr   = (shape_fyle.data(:,1) - min(shape_fyle.data(:,1)))*dtvals;
                    
                    % Compute autocorrelation
                    shape_autocf = compute_autocorr(kappasqout,1);
                    
                    % write autocorrelation output
                    fw_com = fopen(sprintf('../../autocorr/outvals/short_autocorr_gMW_%d_config_%d_sig_%g_eps_%g.dat',...
                        gr_mw,config,sig_val,eps_val),'w');
                    fprintf(fw_com,'%s\t%s\t%s\t%s\n','Step','Time','Autocorr','Normautocorr');
                    fprintf(fw_com,'%g\t%g\t%g\t%g\n',[shape_fyle.data(:,1) time_arr shape_autocf shape_autocf/shape_autocf(1)]');
                    fclose(fw_com);
                    
                end %end of epsilon loop
                
            end %end of sigma loop
            
        end % end of graft MW loop
        
    end %end of bb MW loop
    
end %end of config loop
