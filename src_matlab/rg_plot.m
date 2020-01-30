function rg_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,rg0,rgstdev)

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'k','r',green,'m',brown,'b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};

%% Zero Arrays and Compute
rgavg_vals = zeros(length(num_eps_arr),length(num_sigma_arr));
cntr_arr   = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_rg  = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_rg    = zeros(length(num_eps_arr),length(num_sigma_arr));
rgall_config_vals = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));

for conf_cnt = 1:length(config_arr) %% Config array
    rg_data_all = importdata(sprintf('../../outfiles/config_%d/rgavg_bbMW_%d_gMW_%d_nch_%d.dat', ...
        config_arr(conf_cnt),bb_mw,gr_mw,num_bb_chains));
    
    rgarr_data = rg_data_all.data;
    len_rgdata = length(rgarr_data(:,1));
    
    %map output data to my epsilon/sigma arrays
    for datcnt = 1:len_rgdata
        
        data_epsval = rgarr_data(datcnt,1);
        data_sigval = rgarr_data(datcnt,2);
        
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
        
        rgall_config_vals(myeps_cnt,mysig_cnt,conf_cnt) = rgarr_data(datcnt,3);
        rgavg_vals(myeps_cnt,mysig_cnt) = rgavg_vals(myeps_cnt,mysig_cnt) + rgarr_data(datcnt,3);
        if rgarr_data(datcnt,3) ~= 0 %% add if and only if value is not zero
            cntr_arr(myeps_cnt,mysig_cnt) = cntr_arr(myeps_cnt,mysig_cnt)+1;
        end
        
    end
end

%% Compute Mean and Std Dev
for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        rgavg_vals(eps_cnt,sig_cnt) = rgavg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_avg = zeros(length(config_arr),1);
        
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_avg(conf_cnt,1) = rgall_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        
        % now weed out all the non-zero elements using find command
        [~,~,avgarr] = find(matlab_recheck_avg(:,1));
        
        mean_rg(eps_cnt,sig_cnt)   = mean(avgarr);
        stddev_rg(eps_cnt,sig_cnt) = std(avgarr);
        
    end
end

%% Fig1: Without Error Bar
% Axes Labels and Plot Dimensions
h_avgrg = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$R_g/R_{g0}$','FontSize',20,'Interpreter','Latex')
xlim([0.0 0.02+max(num_sigma_arr)]);
for eps_cnt = 1:length(num_eps_arr)
    plot(num_sigma_arr,rgavg_vals(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{pg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_avgrg,sprintf('./../../all_figures/fig_avgrg_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt;


%% Fig2: Rg With Error Bar
% Axes Labels and Plot Dimensions
h_rgwitherr = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$R_g$','FontSize',20,'Interpreter','Latex')
xlim([0.0 0.02+max(num_sigma_arr)]);
for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,mean_rg(eps_cnt,:),stddev_rg(eps_cnt,:),...
        'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{pg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_rgwitherr,sprintf('./../../all_figures/fig_avgrgwitherr_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt;

%% Fig3: Rg/Rg0 With Error Bar
% Axes Labels and Plot Dimensions
h_rgbyrg0witherr = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$R_g/R_{\rm{g0}}$','FontSize',20,'Interpreter','Latex')
xlim([0.0 0.02+max(num_sigma_arr)]);
for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,mean_rg(eps_cnt,:)/rg0,stddev_rg(eps_cnt,:)./rgstdev,...
        'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{pg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_rgbyrg0witherr,sprintf('./../../all_figures/fig_avgrgbyrg0witherr_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt;

