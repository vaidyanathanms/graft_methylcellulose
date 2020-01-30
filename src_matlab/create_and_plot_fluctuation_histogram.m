function create_and_plot_fluctuation_histogram(shapesrt,nch,lfyle,config,bb_mw,gr_mw,eps_val,sig_val)

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {orange,green,'m','k',brown,'b','r',gold};
lsty = {'-','--',':'};

nbins = 100;
% convert shapesrt into nch arrays of size lfyle/nch*2
shapedata = zeros(floor(lfyle/nch),2,nch);
chcntr    = zeros(nch,1);
for j = 1:lfyle
    chid = shapesrt(j,2);
    chcntr(chid,1) = chcntr(chid,1) + 1;
    shapedata(chcntr(chid,1),1,chid) = shapesrt(j,1);
    shapedata(chcntr(chid,1),2,chid) = shapesrt(j,7);
end

% plot histogram of each chain
hhist = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\kappa^2$','FontSize',16,'Interpreter','Latex')
ylabel('$p(\kappa^2)$','FontSize',16,'Interpreter','Latex')
for i = 1:nch
    histogram(shapedata(:,2,i),nbins,'Normalization','probability');
end

saveas(hhist,sprintf('./../../all_figures/config_%d/fig_hist_bbMW_%d_gMW_%d_nch_%d_eps_%g_sig_%g.png',...
    config,bb_mw,gr_mw,nch,eps_val,sig_val))   


% plot shapefactor for each chain as a function of time
hplot = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('Time','FontSize',16,'Interpreter','Latex')
ylabel('$\kappa^2$','FontSize',16,'Interpreter','Latex')
for i = 1:nch
    plot(shapedata(:,1,i),shapedata(:,2,i),'Color',pclr{i},'LineStyle',lsty{1},'LineWidth',2)
end

saveas(hplot,sprintf('./../../all_figures/config_%d/fig_kappatime_bbMW_%d_gMW_%d_nch_%d_eps_%g_sig_%g.png',...
    config,bb_mw,gr_mw,nch,eps_val,sig_val))   
