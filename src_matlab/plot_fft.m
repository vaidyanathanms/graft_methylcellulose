%% Plot frequency spectra
clc;
clear;
close all;
format long

%% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];
pclr = {orange,green,'m','k',brown,'b','r',gold};
lsty = {'-','--',':'};

%% Input data

sample_fftdata= 1; %resample fft data
num_bb_chains = 1;
polydens      = '0.5';
eps_arr       = {'0.8'};%'0.8';'1.0';'1.2'};
sigma_arr     = {'0.01','0.1','0.2','0.3'};%'0.05','0.1','0.15','0.2','0.25','0.3'};
mw_graft_arr  = [25];
mw_bb_arr     = [1000];
config_arr    = [2];
src_folder   = pwd;

%% FFT Data
%Source: https://blogs.uoregon.edu/seis/wiki/unpacking-the-matlab-fft/
nsamples = 1048*2; %keep it even for nyquist frequency
tminval = 7000;
tsample = 6;
tmaxval = tminval + nsamples*tsample;
tinvals = tminval:tsample:tmaxval;
xdim_data = 2; 
ydim_data = 3; 
ydim_fft  = 4; % IMPORTANT: Most likely this will be ydim_data + 1 because of the way the sample_fft is coded. 

for conf_cnt = 1:length(config_arr)
    config = config_arr(conf_cnt);
    
    for bb_cnt = 1:length(mw_bb_arr) % backbone MW loop
        bb_mw  = mw_bb_arr(bb_cnt);
        
        for gr_cnt = 1:length(mw_graft_arr) %graft MW loop
            gr_mw = mw_graft_arr(gr_cnt);
            
            for eps_cnt = 1:length(eps_arr) %epsilon loop
                eps_val = str2double(eps_arr{eps_cnt});
                
                
                hsq_freq = figure; %with inset
                hold on
                box on
                set(gca,'FontSize',16)
                ax1 = gca; % Store handle to axes 1.
                xlabel(ax1,'Frequency ($\tau^{-1}$)','FontSize',20,'Interpreter','Latex')
                ylabel(ax1,'$\kappa^2 (\omega)$','FontSize',20,'Interpreter','Latex')
                
                %plot inset: https://www.mathworks.com/matlabcentral/answers/60376-how-to-make-an-inset-of-matlab-figure-inside-the-figure?s_tid=gn_loc_drop
                ax2 = axes('Position',[.25 .25 .25 .25]);
                box on;
                    
                
                hsq_time = figure;
                hold on
                box on
                set(gca,'FontSize',20)
                xlabel('Time ($\tau$)','FontSize',20,'Interpreter','Latex')
                ylabel('$\kappa^2$','FontSize',20,'Interpreter','Latex')
                
                for sig_cnt = 1:length(sigma_arr) %sigma loop
                    sig_val = str2double(sigma_arr{sig_cnt});
                    
                    if sample_fftdata %for sampling equispaced intervals
                        fprintf('Sampling time series for config/gr_mw/epsval: %d\t%d\t%g\n',...
                            config,gr_mw,eps_val)
                        
                        fylename = sprintf('../../autocorr/kappasq_time/kappasqtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                            config,gr_mw,sig_val,eps_val); 
                        if exist(fylename,'file') ~= 2
                            fprintf('%s does not exist\n',time_fyle);
                            continue;
                        elseif struct(dir(fylename)).bytes == 0
                            fprintf('Empty file: %s \n',time_fyle);
                            continue;
                        end
                        kappasq = importdata(fylename);    
                        
                        % Sample fft and store sampled data
                        set(0, 'CurrentFigure',hsq_time)
                        ktime = sample_fft(kappasq.data, xdim_data, ydim_data, tinvals, tsample, config, gr_mw, sig_val, eps_val);
                        plot(ktime(:,1),ktime(:,2),'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                        legendtime{sig_cnt} = ['$\Sigma$: ' num2str(sig_val)];
                    end
                    
                    fprintf('Computing and plotting frequency spectra for config/gr_mw/epsval: %d\t%d\t%g\n',...
                        config,gr_mw,eps_val)
                    fylename = sprintf('../../autocorr/fft_samples/fftsampletime_conf_%d_grmw_%d_sig_%g_eps_%g',...
                        config,gr_mw,sig_val,eps_val);
                    if exist(fylename,'file') ~= 2
                        error('%s does not exist\n',time_fyle);
                        continue;
                    elseif struct(dir(fylename)).bytes == 0
                        error('Empty file: %s \n',time_fyle);
                        continue;
                    end
                    
                    sampled_kappasq = importdata(fylename);
                    [fft_data,freq_vec] = compute_fft(sampled_kappasq.data, ydim_fft, tsample, nsamples);
                    
                    index_freqvec = 1:length(freq_vec);                    % Index Vector
                    
                    %write output (half spectra) -- should be symmetric for a
                    %real signal
                    fylename = sprintf('../../autocorr/fft_samples/fftoutput_conf_%d_grmw_%d_sig_%g_eps_%g',...
                        config,gr_mw,sig_val,eps_val);
                    fwrite_fft = fopen(fylename,'w');
                    fprintf(fwrite_fft,'%s\t%s\n','omega','abs(fft)');
                    fprintf(fwrite_fft,'%g\t%g\n',[freq_vec, abs(fft_data(index_freqvec))]'); %two-sided frequency plot
                    fclose(fwrite_fft);
                    
                    
                    % Plot only half the spectra
                    set(0, 'CurrentFigure',hsq_freq)
                    plot(ax1,freq_vec,abs(fft_data(index_freqvec)),'Color',pclr{sig_cnt},...
                        'Marker','x','MarkerSize',8,'LineStyle','None') %double the amplitude
                    legendfreq{sig_cnt} = ['$\Sigma$: ' num2str(sig_val)];
    
                    % Overlay inset
                    plot(ax2, ktime(:,1),ktime(:,2),'Color',pclr{sig_cnt},'LineStyle','-','LineWidth',2)
                    hold on
                    clear kappasq
                    
                end % end of sigma loop
                
                xlim([tminval tmaxval])
                
                set(ax1,'yscale','log');
                set(ax1,'xscale','log');
                legend(ax1,legendfreq,'Interpreter','Latex','FontSize',14,'Location','NorthEast','box','off')
                
                ylabel(ax2,'\kappa^2','FontSize',10)
                xlabel(ax2,'Time (\tau)','FontSize',10)
                
                saveas(hsq_freq,sprintf('./../../all_figures/fig_fftksqtime_conf_%d_gMW_%d_eps_%g.png',...
                    config,gr_mw,eps_val))
                
                set(0, 'CurrentFigure',hsq_time)
                legend(legendtime,'Interpreter','Latex','FontSize',14,'Location','best')
                legend boxoff
                saveas(hsq_time,sprintf('./../../all_figures/fig_rawksqtime_conf_%d_gMW_%d_eps_%g.png',...
                    config,gr_mw,eps_val))
                
            end %end of eps loop
            
        end % end of graft MW loop
        
    end %end of bbMW loop
    
end % end of config

