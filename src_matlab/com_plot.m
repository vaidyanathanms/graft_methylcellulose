function com_plot(config_arr,num_eps_arr,num_sigma_arr,bb_mw,gr_mw,num_bb_chains,dcom_bare,err_dcom_bare)

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {'k',orange,green,'m',brown,'b', gold};
lsty = {'-','--',':'};
msty = {'d','s','o','x'};


%% Plot Mean and Standard Deviations of dCOM over different configurations.

% Zero Arrays and Compute
comavg_vals  = zeros(length(num_eps_arr),length(num_sigma_arr));
cntr_arr     = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_com   = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_norm  = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_com     = zeros(length(num_eps_arr),length(num_sigma_arr));
comall_config_vals = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));

fprintf('Averaging.. \n')

% Main Analysis
for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    fylename = sprintf('../../outfiles/config_%d/distcomavg_bbMW_%d_gMW_%d_nch_%d.dat',...
        config,bb_mw,gr_mw,num_bb_chains);
    
    if exist(fylename,'file') ~= 2
        fprintf('%s \t not found \n',fylename)
        continue;
    end
    
    alldata = importdata(fylename);
    if isstruct(alldata) ~= 1
        fprintf('No data found in %s\n',fylename)
        continue;
    end
    
    comarr_data = alldata.data;
    len_data  = length(comarr_data(:,1));
    
    %map output data to my epsilon/sigma arrays
    for datcnt = 1:len_data
        
        data_sigval = comarr_data(datcnt,1);
        data_epsval = comarr_data(datcnt,2);
        
        find_eps = -1;
        for myeps_cnt = 1:length(num_eps_arr) %epsilon loop
            if data_epsval == num_eps_arr(myeps_cnt)
                find_eps = 1;
                break;
            end
        end
        
        if find_eps == -1
            fprintf('ERROR: Unknown epsilon value in output data: %g for config %g\n',...
                data_epsval,config);
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
            fprintf('ERROR: Unknown sigma value in output data: %g for config %g\n',...
                data_sigval,config);
            break;
        end
        
        comall_config_vals(myeps_cnt,mysig_cnt,conf_cnt) = comarr_data(datcnt,3);
        comavg_vals(myeps_cnt,mysig_cnt) = comavg_vals(myeps_cnt,mysig_cnt) + comarr_data(datcnt,3);
        if comarr_data(datcnt,3) ~= 0 %% add if and only if value is not zero
            cntr_arr(myeps_cnt,mysig_cnt) = cntr_arr(myeps_cnt,mysig_cnt)+1;
        end
        
    end
    
end


%% Compute Mean and Std Error of the Mean

fprintf('Computing mean and standard deviation.. \n');

for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        comavg_vals(eps_cnt,sig_cnt) = comavg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_com_avg = zeros(length(config_arr),1);
        
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_com_avg(conf_cnt,1) = comall_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end
        
        % now weed out all the non-zero elements using find command
        [~,~,avgarr] = find(matlab_recheck_com_avg(:,1));
        
        mean_com(eps_cnt,sig_cnt)   = mean(avgarr);
        stddev_com(eps_cnt,sig_cnt) = std(avgarr)/sqrt(length(avgarr)); %Standard error of the mean is plotted
        
        term1 = (stddev_com(eps_cnt,sig_cnt)/mean_com(eps_cnt,sig_cnt))^2;
        term2 = (err_dcom_bare/dcom_bare)^2;
        stddev_norm(eps_cnt,sig_cnt) = (mean_com(eps_cnt,sig_cnt)/dcom_bare)*sqrt(term1+term2);
        
    end
end

for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        if comavg_vals(eps_cnt,sig_cnt) ~= mean_com(eps_cnt,sig_cnt)
            fprintf('ERROR: @ sig: %g/eps: %g',num_sigma_arr(sig_cnt),num_eps_arr(eps_cnt))
            fprintf('Value from mean_rg: %g\t, rgavg_vals: %g\n', mean_com(eps_cnt,sig_cnt),comavg_vals(eps_cnt,sig_cnt));
            error('ERROR: Some problem in computing averages from two methods: See Sec: Compute Mean and Std Dev\n');
        end
    end
end

%% Plotting Data

fprintf('Plotting average COM over different configurations ..\n');
% Plot Data
h_avgcom = figure;
hold on
box on
set(gca,'FontSize',20)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')


for eps_cnt = 1:length(num_eps_arr)
    plot(num_sigma_arr,comavg_vals(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_avgcom,sprintf('./../../all_figures/fig_avgcom_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))
clear legendinfo eps_cnt sig_cnt alldata;

% Plot with error bar
h_avgerrcom = figure;
hold on
box on
set(gca,'FontSize',20)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,comavg_vals(eps_cnt,:),stddev_com(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
saveas(h_avgerrcom,sprintf('./../../all_figures/fig_avgcomwitherrbar_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt alldata;

%Plot after scaling with bare value
h_normwitherr = figure;
hold on
box on
set(gca,'FontSize',20)
xlabel('$\Sigma$','FontSize',20,'Interpreter','Latex')
ylabel('$d^{*}_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

for eps_cnt = 1:length(num_eps_arr)
    errorbar(num_sigma_arr,comavg_vals(eps_cnt,:)/dcom_bare,stddev_norm(eps_cnt,:),...
        'Color',pclr{eps_cnt},'MarkerSize',12,'MarkerFaceColor',...
        pclr{eps_cnt},'Marker',msty{eps_cnt},'LineStyle',':')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end

set(gca,'yscale','log')
legend(legendinfo,'Interpreter','Latex','FontSize',14,'Location','best')
legend boxoff
saveas(h_normwitherr,sprintf('./../../all_figures/fig_avgcomwitherrbar_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear legendinfo eps_cnt sig_cnt alldata;
%-----------------------------------------------------------------------------------------
%-----------------------------end config averaged plots--------------------



%% Plot dCOM Scaled with Rg

fprintf('Plotting dCOM scaled with Rg ..\n');

% Zero Arrays and Compute Rg Initially
rgavg_vals = zeros(length(num_eps_arr),length(num_sigma_arr));
cntr_arr   = zeros(length(num_eps_arr),length(num_sigma_arr));
stddev_rg  = zeros(length(num_eps_arr),length(num_sigma_arr));
mean_rg    = zeros(length(num_eps_arr),length(num_sigma_arr));
rgall_config_vals = zeros(length(num_eps_arr),length(num_sigma_arr),length(config_arr));

for conf_cnt = 1:length(config_arr) %% Config array
    rg_data_all = importdata(sprintf('../../outfiles/config_%d/rgavg_bbMW_%d_gMW_%d_nch_%d.dat', ...
        config_arr(conf_cnt),bb_mw,gr_mw,num_bb_chains));
    
    if isstruct(rg_data_all) ~= 1
        fprintf('No data for Rg at config %d \n',config_arr(conf_cnt));
        continue;
    end
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
            fprintf('ERROR: Unknown epsilon value in Rg data: %g at config %g\n',...
                data_epsval,config);
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
            fprintf('ERROR: Unknown sigma value in Rg data: %g at config %g\n',...
                data_sigval,config);
            break;
        end
        
        rgall_config_vals(myeps_cnt,mysig_cnt,conf_cnt) = rgarr_data(datcnt,3);
        rgavg_vals(myeps_cnt,mysig_cnt) = rgavg_vals(myeps_cnt,mysig_cnt) + rgarr_data(datcnt,3);
        if rgarr_data(datcnt,3) ~= 0 %% add if and only if value is not zero
            cntr_arr(myeps_cnt,mysig_cnt) = cntr_arr(myeps_cnt,mysig_cnt)+1;
        end
        
    end
end

%% Compute Mean and Std Dev of Rg
for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        rgavg_vals(eps_cnt,sig_cnt) = rgavg_vals(eps_cnt,sig_cnt)/cntr_arr(eps_cnt,sig_cnt);
        matlab_recheck_rg_avg = zeros(length(config_arr),1);
        
        for conf_cnt = 1:length(config_arr)
            matlab_recheck_rg_avg(conf_cnt,1) = rgall_config_vals(eps_cnt,sig_cnt,conf_cnt);
        end

        % now weed out all the non-zero elements using find command
        [~,~,avgarr] = find(matlab_recheck_rg_avg(:,1));
        
        mean_rg(eps_cnt,sig_cnt)    = mean(avgarr);
        stddev_rg(eps_cnt,sig_cnt)  = std(avgarr)/sqrt(length(avgarr));
        
    end
end


%% Plot dCOM as a function of Rg
hrgcom2 = figure;
hold on
box on
set(gca,'FontSize',16)
ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
xlabel('$R_{g}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
for eps_cnt = 1:length(num_eps_arr)
    plot(mean_rg(eps_cnt,:),mean_com(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,...
        'MarkerFaceColor',pclr{eps_cnt},'Marker',msty{3},'LineStyle','None')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')

saveas(hrgcom2,sprintf('./../../all_figures/fig_scalingdandrgdata_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))
clear legendinfo

xstd = stddev_rg;
ystd = stddev_com;

%% Plot dCOM as a function of Rg with errorbar
hrgcom2 = figure;
hold on
box on
set(gca,'FontSize',16)

ylabel('$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
xlabel('$R_{g}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
for eps_cnt = 1:length(num_eps_arr)
    errorbar(mean_rg(eps_cnt,:),mean_com(eps_cnt,:),ystd(eps_cnt,:),ystd(eps_cnt,:),...
        xstd(eps_cnt,:),xstd(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,...
        'MarkerFaceColor',pclr{eps_cnt},'Marker',msty{3},'LineStyle','None')
    legendinfo{eps_cnt} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
end
legend(legendinfo,'Interpreter','Latex','FontSize',20,'Location','best')
legend boxoff
saveas(hrgcom2,sprintf('./../../all_figures/fig_scalingdandrgdatawitherr_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear xdata ymeandata ystddata rgoutdata legendinfo


%% Plot dCOM with and without errorbar (inset)
hinset = figure;
hold on
box on
set(gca,'FontSize',16)
ax1 = gca; % Store handle to axes 1.
ylabel(ax1,'$d_{\rm{COM}}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')
xlabel(ax1,'$R_{g}$ ($\sigma$)','FontSize',20,'Interpreter','Latex')

%plot inset: https://www.mathworks.com/matlabcentral/answers/60376-how-to-make-an-inset-of-matlab-figure-inside-the-figure?s_tid=gn_loc_drop
ax2 = axes('Position',[.7 .7 .25 .25]);
box on;

% copy all data into a single array for line fit

fulldata = zeros(length(num_eps_arr)*length(num_sigma_arr),2);
cntr = 0;
for eps_cnt = 1:length(num_eps_arr)
    for sig_cnt = 1:length(num_sigma_arr)
        cntr = cntr + 1;
        fulldata(cntr,1) = mean_rg(eps_cnt,sig_cnt);
        fulldata(cntr,2) = mean_com(eps_cnt,sig_cnt);
    end
end
   
lin_coeffs = polyfit(fulldata(:,1),fulldata(:,2),1);
% Get fitted values
fittedX = linspace(0, 1.5*max(fulldata(:,1)), 200);
fittedY = polyval(lin_coeffs, fittedX);

plot(ax1,fittedX, fittedY, 'b--', 'LineWidth', 2);
legendinfo{1} = ['Fitted Line'];
hold on

for eps_cnt = 1:length(num_eps_arr)
    errorbar(ax1,mean_rg(eps_cnt,:),mean_com(eps_cnt,:),ystd(eps_cnt,:),ystd(eps_cnt,:),...
        xstd(eps_cnt,:),xstd(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',12,...
        'MarkerFaceColor',pclr{eps_cnt},'Marker',msty{3},'LineStyle','None')
    legendinfo{eps_cnt+1} = ['$\epsilon_{gg}$: ' num2str(num_eps_arr(eps_cnt))];
    % Overlay inset
    plot(ax2,mean_rg(eps_cnt,:),mean_com(eps_cnt,:),'Color',pclr{eps_cnt},'MarkerSize',3,...
        'MarkerFaceColor',pclr{eps_cnt},'Marker',msty{3},'LineStyle','None')
    hold on
end



xlim(ax1,[0 1.5*max(max(mean_rg+xstd))]);
ylim(ax1,[0 1.5*max(max(mean_com+ystd))]);
xlabel(ax2,'$R_{g}$ ($\sigma$)','FontSize',10,'Interpreter','Latex')
ylabel(ax2,'$d_{\rm{COM}}$ ($\sigma$)','FontSize',10,'Interpreter','Latex')

legend(ax1,legendinfo,'Interpreter','Latex','FontSize',14,'Location','NorthWest','box','off')
legend boxoff
legend(ax2,'None');
legend boxoff

saveas(hinset,sprintf('./../../all_figures/fig_scalingdandrgdatainset_bbMW_%d_gMW_%d_nch_%d.png',...
    bb_mw,gr_mw,num_bb_chains))

clear xdata ymeandata ystddata rgoutdata legendinfo