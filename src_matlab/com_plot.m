function com_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains)

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'k','r',green,'m',brown,'b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};


%% Plot Mean and Standard Deviations of dCOM over different configurations.

% Zero Arrays and Compute
comavg_vals  = zeros(length(num_eps_arr),length(num_sigma_arr));
cntr_arr     = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_com   = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_com     = zeros(length(num_eps_arr),length(num_sigma_arr));
comall_config_vals = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));

% Axes Labels and Plot Dimensions
h_avgcom = figure;
hold on
box on
set(gca,'FontSize',20)
xlabel('\Sigma','FontSize',20,'Interpreter','Latex')
ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

% Main Analysis
for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    alldata = importdata(sprintf('../../outfiles/config_%d/distcom_bbMW_%d_gMW_%d_nch_%d.dat',...
        config,bb_mw,gr_mw,num_bb_chains));
    
    comarr_data = alldata.data;
    len_data  = length(comarr_data(:,1));
    
    %map output data to my epsilon/sigma arrays
    for datcnt = 1:len_rgdata
        
        data_epsval = comarr_data(datcnt,1);
        data_sigval = comarr_data(datcnt,2);
        
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
        
        comall_config_vals(myeps_cnt,mysig_cnt,conf_cnt) = comarr_data(datcnt,4);
        comavg_vals(myeps_cnt,mysig_cnt) = comavg_vals(myeps_cnt,mysig_cnt) + comarr_data(datcnt,3);
        if comarr_data(datcnt,3) ~= 0 %% add if and only if value is not zero
            cntr_arr(myeps_cnt,mysig_cnt) = cntr_arr(myeps_cnt,mysig_cnt)+1;
        end
        
    end
    
end


% Compute Mean and Std Dev
for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        comavg_vals(eps_cnt,sig_cnt) = comavg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_avg = zeros(length(config_arr),1);
        
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_avg(conf_cnt,1) = rgall_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        
        mean_com(eps_cnt,sig_cnt)   = mean(matlab_recheck_avg(:,1));
        stddev_com(eps_cnt,sig_cnt) = std(matlab_recheck_avg(:,1));
        
    end
end

% Plot Data
for eps_cnt = 1:length(num_eps_arr)
    plot(num_sigma_arr,comavg_vals(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{pg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_avgcom,sprintf('./../../all_figures/fig_avgcom_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt alldata;

%-----------------------------------------------------------------------------------------
%-----------------------------end config averaged plots--------------------



%% Plot each configuration independently

% Axes Labels and Plot Dimensions
h2com = figure;
hold on
box on
set(gca,'FontSize',20)
xlabel('\Sigma','FontSize',20,'Interpreter','Latex')
ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

for conf_cnt = 1:length(config_arr)
    
    config = config_arr(conf_cnt);
    alldata = importdata(sprintf('../../outfiles/config_%d/distcom_bbMW_%d_gMW_%d_nch_%d.dat',...
        config,bb_mw,gr_mw,num_bb_chains));
    comarr_data = alldata.data;
    
    for eps_cnt = 1:length(num_eps_arr) %epsilon loop
        
        xdata     = zeros(length(num_sigma_arr),1);
        ymeandata = zeros(length(num_sigma_arr),1);
        ystddata  = zeros(length(num_sigma_arr),1);
        sig_cnt = 0;
        eps_val = num_eps_arr(eps_cnt);
        for fcnt = 1:len_data
            if comarr_data(fcnt,2) == eps_val
                sig_cnt = sig_cnt + 1;
                xdata(sig_cnt,1)     = comarr_data(fcnt,1);
                ymeandata(sig_cnt,1) = comarr_data(fcnt,4);
                ystddata(sig_cnt,1)  = comarr_data(fcnt,5);
            end
        end
        plot(xdata,ymeandata,'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
            pclr{eps_cnt},'Marker',msty{1},'LineStyle',':')
        legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
    end
    legendinfo{eps_cnt} = ['$\epsilon_{pg}$: ' num2str(num_eps_arr(eps_cnt))];
    saveas(h2com,sprintf('./../../all_figures/config_%d/fig_comdata_bbMW_%d_gMW_%d_nch_%d_%s.png',...
        config,bb_mw,gr_mw,num_bb_chains,num2str(num_eps_arr(eps_cnt))));
end 

clear legendinfo eps_cnt sig_cnt

%--------------------------------------------------------------------------
%---------------------------End of independent configs plot----------------

%% Plot dCOM Scaled with Rg

hrgcom2 = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$d_{\rm{COM}$','FontSize',20,'Interpreter','Latex')
ylabel('R_{g}','FontSize',20,'Interpreter','Latex')

% Zero Arrays and Compute Rg Initially
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

% Compute Mean and Std Dev of Rg
for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        rgavg_vals(eps_cnt,sig_cnt) = rgavg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_avg = zeros(length(config_arr),1);
        
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_avg(conf_cnt,1) = rgall_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        
        mean_rg(eps_cnt,sig_cnt)   = mean(matlab_recheck_avg(:,1));
        stddev_rg(eps_cnt,sig_cnt) = std(matlab_recheck_avg(:,1));
        
    end
end

for eps_cnt = 1:length(num_eps_arr)
    plot(mean_com(eps_cnt,:),mean_rg(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,...
        'MarkerFaceColor',pclr{eps_cnt},'Marker',msty{3},'LineStyle',':')
end

% Plot dCOM as a function of Rg

saveas(hrgcom2,sprintf('./../../all_figures/fig_scalingdandrgdata_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))
clear xdata ymeandata ystddata rgoutdata legendinfo