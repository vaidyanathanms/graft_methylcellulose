%% To analyze static properties

clc;
clear;
close all;
format long;

%% Input Flags

rgflag   = [0,0,1,2];%[write,plot,XY-analysis_column]
persflag = [0,0,1,2];%[write,plot,XY-analysis_column]
comflag  = [0,1,1,5];%[write,plot,XY-analysis_column]
shapeflag = [0,0];%[write,plot]

%% Input Data

num_bb_chains = 2;
polydens      = '0.1';
eps_arr       = {'0.8','1.0','1.2'};
sigma_arr     = {'0.01','0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [1,2,4,5]; 

%% Bare chain(s) results
rg_bare       = 9.78; %mean Rg of bare systems - no graft
err_rg_bare   = 0.25; %Reporting SEM; std dev = 0.60 for bare systems
dcom_bare     = 3.61193;
err_dcom_bare = 0.59; %Reporting SEM; std dev = 1.69109 for bare systems

%% Cut-off data
cutoff       = 0.1; %cutoff for persistence length
tcut_frac    = 0.5; %cutoff time fraction to start calculating dcomavg

mincut = 3*10^7; %cutoff timestep beyond which avg values are computed
maxcut = 7*10^7; %maximum timestep over which averages are computed

%% Main analysis - Write consolidated data

for conf_cnt = 1:length(config_arr) %config loop (replicate trials)
    
    config = config_arr(conf_cnt);
    fprintf('Configuration under analysis: %d\n', config);
    
    for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
    
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
        
            %Initiate output files and get file ID
            if rgflag(1,1)
                fw_rg = fopen(sprintf('../../outfiles/config_%d/rgavg_bbMW_%d_gMW_%d_nch_%d.dat',...
                    config,bb_mw,gr_mw,num_bb_chains),'w');
                fprintf(fw_rg,'%s\t%s\t%s\t%s\t%s\n','epsilon','sigma','Avg_rg','cutoff_min_step','cutoff_max_step');
            end
            
            if persflag(1,1)
                fw_per = fopen(sprintf('../../outfiles/config_%d,pers_bbMW_%d_gMW_%d_nch_%d.dat',...
                    config,bb_mw,gr_mw,num_bb_chains),'w');
                fprintf(fw_per,'%s\n','Fit fun y = cf1*exp(-cf2*x)');
                fprintf(fw_per,'%s\t%s\t%s\t%s\t%s\t%s\n','sigma','epsilon',...
                    'cf1','lp=1/cf2','npoints','fitquality');
            end
            
            if comflag(1,1)
                fw_com = fopen(sprintf('../../outfiles/config_%d/distcomavg_bbMW_%d_gMW_%d_nch_%d.dat',...
                    config,bb_mw,gr_mw,num_bb_chains),'w');
                fprintf(fw_com,'%s\t%s\t%s\t%s\t%s\t%s\n','sigma','epsilon','mean_COM','stddev_COM','cutoff_min_step','cutoff_max_step');
            end
            
            if shapeflag(1,1)
                fw_shape = fopen(sprintf('../../outfiles/config_%d/shapefac_bbMW_%d_gMW_%d_nch_%d.dat',...
                    config,bb_mw,gr_mw,num_bb_chains),'w');
                
                fprintf(fw_shape,'%s\t %s%d%s\t','sigma','epsilon (Cutoff_Time: ',mincut,')');
                for chcnt = 1:num_bb_chains
                    rg_str   = sprintf('%s%d','avg_rg: ',chcnt);
                    sph_str  = sprintf('%s%d','avg_sphericity: ',chcnt);
                    cyl_str  = sprintf('%s%d','avg_cylindricity: ',chcnt);
                    shp_str  = sprintf('%s%d','avg_shape: ',chcnt);
                    prol_str = sprintf('%s%d','avg_prolateness: ',chcnt);
                    fprintf(fw_shape, '%s\t%s\t%s\t%s\t%s\t',rg_str,sph_str,cyl_str,...
                        shp_str,prol_str);
                end
                fprintf(fw_shape,'\n');
                clear cut_str sph_str cyl_str shp_str
            end
            
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
                    
                    if rgflag(1,1) %begin rg calculation
                        fprintf('Analyzing Rg values \n');
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
                            
                            [mincutoff,maxcutoff] = find_general_cutoff_minmax(sortrg(:,1),mincut,maxcut);
                             
                            rg_avg = mean(sortrg(mincutoff:maxcutoff,2));
                            fprintf(fw_rg,'%g\t%g\t%g\t%g\t%g\n',eps_val,sig_val,rg_avg,sortrg(mincutoff,1),sortrg(maxcutoff,1));
                            clear xcol ycol mincutoff maxcutoff
                        end
                        clear rg_fylename
                    end % End of rg calculation
                    
                    if persflag(1,1) %begin persistence length calculation
                        fprintf('Analyzing persistence length values \n');
                        per_prefix = sprintf('mainpersistautocf_%s_%s_%d.dat',...
                            eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                        
                        per_fylename = strcat(dirname,'/',per_prefix);
                        
                        if exist(per_fylename,'file') ~= 2
                            fprintf('%s does not exist\n',per_fylename);
                            continue;
                        elseif struct(dir(per_fylename)).bytes == 0
                            fprintf('Empty file: %s \n',per_fylename);
                            continue;
                        else
                            persdata = importdata(per_fylename);
                            %find cut-off for fit
                            xcol = persflag(1,3); ycol = persflag(1,4);
                            allydata  = persdata.data(:,ycol);
                            cut_index = find_cutindex(allydata,cutoff);
                            %fit data
                            xdata = persdata.data(1:cut_index,xcol);
                            ydata = persdata.data(1:cut_index,ycol);
                            guessvals = [1,1/cut_index];
                            ftype = fittype('a*exp(-b*x)');
                            [coeff,gof] = fit(xdata,ydata,ftype,'Start',guessvals);
                            fprintf(fw_per,'%g\t%g\t%g\t%g\t%g\t%g\n',sig_val,...
                                eps_val,coeff.a,1/coeff.b,cut_index,gof.rsquare);
                            clear xcol ycol
                        end
                        clear per_fylename
                    end % End of persistence length calculation
                    
                    if comflag(1,1) %begin com calculation
                        fprintf('Analyzing com distance values \n');
                        com_prefix = sprintf('composallappend_%s_%s_%d.dat',...
                            eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw);
                        
                        com_fylename = strcat(dirname,'/',com_prefix);
                        
                        if exist(com_fylename,'file') ~= 2
                            fprintf('%s does not exist\n',com_fylename);
                            continue;
                        elseif struct(dir(com_fylename)).bytes == 0
                            fprintf('Empty file: %s \n',com_fylename);
                            continue;
                        else
                            xcol = comflag(1,3); ycol = comflag(1,4);
                            %sort data so that time jumbling is taken care of
                            com_allvals = importdata(com_fylename);
                            com_allvals_unsrt = [com_allvals.data(:,xcol) com_allvals.data(:,ycol)];
                            com_allvals_srted = sortrows(com_allvals_unsrt,xcol);
                            
                            lfyle = length(com_allvals_srted(:,1));
                            [tcut_min, tcut_max] = find_general_cutoff_minmax(com_allvals_srted(:,1),mincut,maxcut);
                            
                            dcom_mean = mean(com_allvals_srted(tcut_min:lfyle,2));
                            dcom_std  = std(com_allvals_srted(tcut_min:lfyle,2));
                            fprintf(fw_com,'%g\t%g\t%g\t%g\t%g\t%g\n',sig_val,eps_val,dcom_mean,dcom_std,...
                                com_allvals_srted(tcut_min,1),com_allvals_srted(tcut_max,1));
                            clear xcol ycol tcut_index com_allvals_srted com_allvals_unsrt com_allvals tcut_min tcut_max
                        end
                        clear com_fylename
                    end %End com calculation
                    
                    if shapeflag(1,1) %begin shapefactor calculation
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
                            shapesrt = sortrows(shapeunsrt,1);
                            
                            lenkappa   = length(shapesrt(:,1));
                            tcut_min = find_general_cutoff_index(shapesrt(:,1),mincut);
                            avg_rg     = zeros(num_bb_chains,1);
                            avg_sphere = zeros(num_bb_chains,1);
                            avg_cylndr = zeros(num_bb_chains,1);
                            avg_shape  = zeros(num_bb_chains,1);
                            avg_prolate = zeros(num_bb_chains,1);
                            ncntr = zeros(num_bb_chains,1);
                            
                            for cntr = tcut_min:lfyle
                                chid = shapesrt(cntr,2);
                                avg_rg(chid,1) = avg_rg(chid,1) + shapesrt(cntr,4);
                                avg_sphere(chid,1) = avg_sphere(chid,1) + shapesrt(cntr,5);
                                avg_cylndr(chid,1) = avg_cylndr(chid,1) + shapesrt(cntr,6);
                                avg_shape(chid,1) = avg_shape(chid,1) + shapesrt(cntr,7);
                                avg_prolate(chid,1) = avg_prolate(chid,1) + shapesrt(cntr,8);
                                ncntr(chid,1) = ncntr(chid,1) + 1;
                            end
            
                            for chainID = 1:num_bb_chains
                                rg_mean    = avg_rg(chainID,1)/ncntr(chainID,1);
                                sph_mean   = avg_sphere(chainID,1)/ncntr(chainID,1);
                                cyl_mean   = avg_cylndr(chainID,1)/ncntr(chainID,1);
                                shape_mean = avg_shape(chainID,1)/ncntr(chainID,1);
                                prol_mean  = avg_prolate(chainID,1)/ncntr(chainID,1);                               
                                fprintf(fw_shape,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t',...
                                    sig_val,eps_val,rg_mean,sph_mean,cyl_mean,shape_mean,prol_mean);
                            end
                            fprintf(fw_shape,'\n');
                            clear xcol ycol tcut_index
                        end
                        clear shape_fylename
                        clear avg_sphere avg_cylndr avg_shape
                    end %end shapefac calculation
                    
                end %end epsilon loop
                
            end %end sigma loop
            
            %Close all files
            if rgflag(1,1)
                fclose(fw_rg);
            end
            if persflag(1,1)
                fclose(fw_per);
            end
            if comflag(1,1)
                fclose(fw_com);
            end
            if shapeflag(1,1)
                fclose(fw_shape);
            end
            
        end %end graftMW
        
    end %end bbMW

end %end config

%% Main analysis - Plot all data

num_sigma_arr = str2double(sigma_arr);
num_eps_arr   = str2double(eps_arr);

for bb_cnt = 1:length(mw_bb_arr) %backbone MW loop
    bb_mw  = mw_bb_arr(bb_cnt);
    
    for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
        gr_mw = mw_graft_arr(gr_cnt);
        
        %% Begin Rg plots
        if rgflag(1,2)
    
            %assumes that the epsilon column is second and sigma column is
            %third
            
            rg_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,rg_bare,err_rg_bare)
            
        end
            
        %% Begin COM plots
        if comflag(1,2) 
            
            com_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,dcom_bare,err_dcom_bare)            
            
        end %End of com plots
        
        %% Begin shape factor plots
        % Assumes the following file format: sigma, eps, avg_rg, avg_sphericity, avg_cylind, avg_shape, avg_prolateness
        if shapeflag(1,2)
                   
            shape_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,plotflag)
                
        end %End of shape plots
            
    end %end gr MW loop
    
end % end bb MW loop


