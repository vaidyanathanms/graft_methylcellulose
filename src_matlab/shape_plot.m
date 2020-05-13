function shape_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains)
% Assumes the following file format: sigma, eps, avg_rg, avg_sphericity, avg_cylind, avg_shape, avg_prolateness

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'k',orange,green,'m',brown,'b','r',gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Zero Arrays and Compute
rg_avg_vals       = zeros(length(num_eps_arr),length(num_sigma_arr));
acyl_avg_vals     = zeros(length(num_eps_arr),length(num_sigma_arr));
asph_avg_vals     = zeros(length(num_eps_arr),length(num_sigma_arr));
prol_avg_vals     = zeros(length(num_eps_arr),length(num_sigma_arr));
shape_avg_vals    = zeros(length(num_eps_arr),length(num_sigma_arr));
cntr_arr          = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_shape      = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_shape        = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_prol       = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_prol         = zeros(length(num_eps_arr),length(num_sigma_arr));
shape_config_vals = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));
prol_config_vals  = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));

for conf_cnt = 1:length(config_arr) %% Config array
    dirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
        config,bb_mw,gr_mw,polydens,num_bb_chains);
    if ~exist(dirname,'dir')
        fprintf('%s does not exist\n',dirname);
        continue
    end
    
    shape_file_data = importdata(sprintf('../../outfiles/config_%d/shapefac_bbMW_%d_gMW_%d_nch_%d.dat',...
        config/bb_mw,gr_mw,num_bb_chains));
    shape_data = shape_file_data.data;
    len_shapedata = length(shape_data);
    
    %map output data to my epsilon/sigma arrays
    for datcnt = 1:len_shapedata
        
        data_sigval = shape_data(datcnt,1);
        data_epsval = shape_data(datcnt,2);
        
        find_eps = -1;
        for myeps_cnt = 1:length(num_eps_arr) %epsilon loop
            if data_epsval == num_eps_arr(myeps_cnt)
                find_eps = 1;
                break;
            end
        end
        
        if find_eps == -1
            fprintf('Unknown epsilon value in output data: %g\n',data_epsval);
            break;
        end
        
        find_sig = -1;
        for mysig_cnt = 1:length(num_sigma_arr) %sigma loop
            if data_sigval == num_sigma_arr(mysig_cnt)
                find_sig = 1;
                break;
            end
        end
        
        if find_sig == -1
            fprintf('Unknown epsilon value in output data: %g\n',data_sigval);
            break;
        end
        
        
        rg_avg_vals(myeps_cnt,mysig_cnt)    = rg_avg_vals(myeps_cnt,mysig_cnt)    + shape_data(datcnt,3);
        acyl_avg_vals(myeps_cnt,mysig_cnt)  = acyl_avg_vals(myeps_cnt,mysig_cnt)  + shape_data(datcnt,4);
        asph_avg_vals(myeps_cnt,mysig_cnt)  = asph_avg_vals(myeps_cnt,mysig_cnt)  + shape_data(datcnt,4);
        shape_avg_vals(myeps_cnt,mysig_cnt) = shape_avg_vals(myeps_cnt,mysig_cnt) + shape_data(datcnt,6);
        prol_avg_vals(myeps_cnt,mysig_cnt)  = prol_avg_vals(myeps_cnt,mysig_cnt)  + shape_data(datcnt,7);
        
        shape_config_vals(myeps_cnt,mysig_cnt,conf_cnt) = shape_data(datcnt,6);
        prol_config_vals(myeps_cnt,mysig_cnt,conf_cnt)  = shape_data(datcnt,7);
        
        
        if shape_data(datcnt,3) ~= 0 %% add if and only if value is not zero -- this has to be always true by defn of code
            cntr_arr(myeps_cnt,mysig_cnt) = cntr_arr(myeps_cnt,mysig_cnt)+1;
        end
        
    end
end

%% Compute mean and stddev 

for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        % avg and std dev of shape factors
        shape_avg_vals(eps_cnt,sig_cnt) = shape_avg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_avg = zeros(length(config_arr),1);
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_avg(conf_cnt,1) = shape_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        % now weed out all the non-zero elements using find command
        [~,~,avgarr] = find(matlab_recheck_avg(:,1));
        mean_shape(eps_cnt,sig_cnt)   = mean(avgarr);
        stddev_shape(eps_cnt,sig_cnt) = std(avgarr);
        
        %avg and std dev of prolateness
        prol_avg_vals(eps_cnt,sig_cnt) = prol_avg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_avg = zeros(length(config_arr),1);
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_avg(conf_cnt,1) = prol_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        % now weed out all the non-zero elements using find command
        [~,~,avgarr] = find(matlab_recheck_avg(:,1));
        mean_prol(eps_cnt,sig_cnt)   = mean(avgarr);
        stddev_prol(eps_cnt,sig_cnt) = std(avgarr);
        
        
    end
end

%% Plot average shape factor with errorbar
% Axes Labels and Plot Dimensions
h_kaperr = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$R_g$','FontSize',20,'Interpreter','Latex')
xlim([0.0 0.02+max(num_sigma_arr)]);
for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,mean_shape(eps_cnt,:),stddev_shape(eps_cnt,:),...
        'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_kaperr,sprintf('./../../all_figures/fig_avgshapewitherr_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt;

%% Plot average prolateness with errorbar
% Axes Labels and Plot Dimensions
h_prolerr = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$R_g$','FontSize',20,'Interpreter','Latex')
xlim([0.0 0.02+max(num_sigma_arr)]);
for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,mean_prol(eps_cnt,:),stddev_prol(eps_cnt,:),...
        'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_prolerr,sprintf('./../../all_figures/fig_avgprolwitherr_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt;


