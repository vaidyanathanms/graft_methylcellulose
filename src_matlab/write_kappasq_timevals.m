function write_kappasq_timevals(kappa_time, kappasqvals, sim_timedata, config,gr_mw,sig_val,eps_val)

fprintf('Writing Kappasq-Time Data ..\n');
fout_kappa  = fopen(sprintf('../../autocorr/kappasq_time/kappasqtime_conf_%d_grmw_%d_sig_%g_eps_%g',...
    config,gr_mw,sig_val,eps_val),'w');

fprintf(fout_kappa,'%s\t%s %s\n','Step','Time','Kappasq');

fout_repeat = fopen(sprintf('../../autocorr/kappasq_time/repeattime_conf_%d_grmw_%d_sig_%g_eps_%g',...
    config,gr_mw,sig_val,eps_val),'w');

fprintf(fout_repeat,'%s\n','Repeated Step');

%% Check whether there are duplicate times
% sort sim_timedata -- if tend has a value greater than the
% succeeding tstart value, multiple values are likely to be written.
lensimtime = length(sim_timedata(:,1));
flags  = zeros(lensimtime,1);
for i = 2:lensimtime
    if sim_timedata(i,1) < sim_timedata(i-1,2)
        flags(i,1) = 1;
        fprintf('WARNING: Check correctness of kappasq data \n');
        fprintf('Found repeat for timesequence between %g\t and %g\t in config/gr_mw/sig_val/eps_val: %d\t%d\t%g\t%g\n',...
            sim_timedata(i,1),sim_timedata(i-1,2),config,gr_mw,sig_val,eps_val);
        fprintf(fout_repeat,'%s\n','WARNING: Check correctness of kappasq data \n');
        fprintf(fout_repeat,'Found repeat for timesequence between %g\t and %g\t in config/gr_mw/sig_val/eps_val: %d\t%d\t%g\t%g\n',...
            sim_timedata(i,1),sim_timedata(i-1,2),config,gr_mw,sig_val,eps_val);
    end
end

%% Compute the corresponding timevalues
telapse = 0; tinit = kappa_time(1,1);
fprintf(fout_kappa,'%g\t%g\t%g\n',kappa_time(1,1),telapse,kappasqvals(1,1));
for j = 2:length(kappa_time)
    % check whether the next time point is same as this - if so skip this
    if j ~= length(kappa_time)
        if kappa_time(j,1) == kappa_time(j+1,1)
            fprintf(fout_repeat,'%g\n',kappa_time(j,1));
            fprintf('Repeat time at index/val(index)/val(index+1): %g\t%g\t%g\n',j, kappa_time(j,1),kappa_time(j+1,1));
            continue;
        end
    end
    % compare with simulation details and write time data
    for k = 1:lensimtime
        if kappa_time(j,1) >= sim_timedata(k,1) && kappa_time(j,1) < sim_timedata(k,2)
            telapse = telapse + (kappa_time(j,1)-tinit)*sim_timedata(k,3); 
            tinit   = kappa_time(j,1);
            fprintf(fout_kappa,'%g\t%g\t%g\n',kappa_time(j,1),telapse,kappasqvals(j,1));
        end
    end
end

fclose(fout_kappa);
fclose(fout_repeat);
            
            

