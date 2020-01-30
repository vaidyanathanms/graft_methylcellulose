function shape_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,plotflags)

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
    dirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
        config,bb_mw,gr_mw,polydens,num_bb_chains);
    if ~exist(dirname,'dir')
        fprintf('%s does not exist\n',dirname);
        continue
    end
    
    shape_file_data = importdata(sprintf('../../outfiles/config_%d/shapefac_bbMW_%d_gMW_%d_nch_%d.dat',...
        config/bb_mw,gr_mw,num_bb_chains));
    shape_data = shape_file_data.data;
    
    
    rgarr_data = rg_data_all.data;
    len_rgdata = length(rgarr_data(:,1));
    
    
    hsph = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('Graft Density','FontSize',16,'Interpreter','Latex')
    ylabel('Sphericity','FontSize',16,'Interpreter','Latex')
    for chcnt = 1:num_bb_chains
        plot(x_pl_data(:,1),sph_pl_data(:,chcnt),'Color',pclr{chcnt},...
            'MarkerSize',8,'MarkerFaceColor',...
            pclr{chcnt},'Marker',msty{1},'LineStyle',':')
        legendinfo{chcnt} = ['Chain ID: ' num2str(chcnt)];
    end
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
    saveas(hsph,sprintf('../../all_figures/config_%d/fig_sphdata_bbMW_%d_gMW_%d_nch_%d_%s.png',...
        config,bb_mw,gr_mw,num_bb_chains,eps_arr{eps_cnt}));
    
    hcyl = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('Graft Density','FontSize',16,'Interpreter','Latex')
    ylabel('Cylindricity','FontSize',16,'Interpreter','Latex')
    for chcnt = 1:num_bb_chains
        
        plot(x_pl_data(:,1),cyl_pl_data(:,chcnt),'Color',pclr{chcnt},...
            'MarkerSize',8,'MarkerFaceColor',...
            pclr{chcnt},'Marker',msty{1},'LineStyle',':')
        legendinfo{chcnt} = ['Chain ID: ' num2str(chcnt)];
    end
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
    saveas(hcyl,sprintf('../../all_figures/config_%d/fig_cyldata_bbMW_%d_gMW_%d_nch_%d_%s.png',...
        config,bb_mw,gr_mw,num_bb_chains,eps_arr{eps_cnt}));
    
    hkap = figure;
    hold on
    box on
    set(gca,'FontSize',16)
    xlabel('Graft Density','FontSize',16,'Interpreter','Latex')
    ylabel('$\langle \kappa^2 \rangle$','FontSize',16,'Interpreter','Latex')
    for chcnt = 1:num_bb_chains
        plot(x_pl_data(:,1),kap_pl_data(:,chcnt),'Color',pclr{chcnt},...
            'MarkerSize',8,'MarkerFaceColor',...
            pclr{chcnt},'Marker',msty{1},'LineStyle',':')
        legendinfo{chcnt} = ['Chain ID: ' num2str(chcnt)];
    end
    legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
    saveas(hkap,sprintf('../../all_figures/config_%d/fig_kapdata_bbMW_%d_gMW_%d_nch_%d_%s.png',...
        config,bb_mw,gr_mw,num_bb_chains,eps_arr{eps_cnt}));
    
    
end

%% Plot individual shapefactors as a function of time
dirname = sprintf('../../sim_results/config_%d/out_bbMW_%d_ngMW_%d_rho_%s_nch_%d',...
    config,bb_mw,gr_mw,polydens,num_bb_chains);
if ~exist(dirname,'dir')
    fprintf('%s does not exist\n',dirname);
    continue
end
for eps_cnt = 1:length(eps_arr) %epsilon loop
    for sig_cnt = 1:2:length(sigma_arr) %sigma loop
        hkaptime = figure;
        hold on
        box on
        set(gca,'FontSize',16)
        xlabel('Timestep','FontSize',16,'Interpreter','Latex')
        ylabel('$\langle \kappa^2 \rangle$','FontSize',16,'Interpreter','Latex')
        for chcnt = 1:num_bb_chains %chain loop
            
            ind_shape_prefix = sprintf('shapeallappend_%s_%s_%d_%d.dat',...
                eps_arr{eps_cnt},sigma_arr{sig_cnt},bb_mw,chcnt);
            ind_shape_fylename = strcat(dirname,'/',ind_shape_prefix);
            ind_shape_alldata = importdata(ind_shape_fylename);
            ind_shape_data = ind_shape_alldata.data;
            xplot_data = ind_shape_data(:,1);
            kap_plot_data = ind_shape_data(:,7);
            plot(xplot_data,kap_plot_data,'Color',pclr{chcnt},'LineStyle','--','LineWidth',2)
            legendinfo{chcnt} = ['Chain ID: ' num2str(chcnt)];
        end
        legend(legendinfo,'Interpreter','Latex','FontSize',16,'Location','best')
        saveas(hkaptime,sprintf('./../../all_figures/config_%d/fig_kapdata_bbMW_%d_gMW_%d_nch_%d_sig_%s_eps_%s.png',...
            config,bb_mw,gr_mw,num_bb_chains,sigma_arr{sig_cnt},eps_arr{eps_cnt}));
    end
end
