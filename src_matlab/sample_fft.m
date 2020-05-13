function ksq_time = sample_fft(indata, xindim, yindim, tinarr, tsample, config, gr_mw, sig_val, eps_val)

%% Sample uniform data from input
tminval = min(tinarr); tmaxval = max(tinarr); ksq_time = zeros(100,2);

% store sampled values in a separate file
fylename = sprintf('../../autocorr/fft_samples/fftsampletime_conf_%d_grmw_%d_sig_%g_eps_%g',...
    config,gr_mw,sig_val,eps_val);
fwrite_time = fopen(fylename,'w');
fprintf(fwrite_time,'%s\t %s\t %s\t %s %g\n','Step (Unshifted)','Unshifted time','Shifted time','Kappasq with sample time: ',tsample);

% Find minimum value and then start from that point while subtracting the
% min time value (DO NOT subtract the step value)
for tfind = 1:length(indata(:,xindim))
    if indata(tfind,xindim) >= tminval
        tmin_index = tfind;
        time_minval = indata(tfind,xindim);
        fprintf(fwrite_time,'%g\t%g\t%g\t%g\n',indata(tfind,1),...
            indata(tfind,xindim),indata(tfind,xindim)-time_minval,indata(tfind,yindim));
        break;
    end
end

curr_time = time_minval;
ksq_time(1,1) = curr_time;
ksq_time(1,2) = indata(tfind,yindim);
arr_index     = 2;
curr_time     = indata(tfind,xindim) + tsample;

for tfind = tmin_index+1:length(indata(:,xindim))
    if indata(tfind,xindim) <= tmaxval %% Sample only at equispaced intervals
        if curr_time == indata(tfind,xindim) 
            fprintf(fwrite_time,'%g\t%g\t%g\t%g\n',indata(tfind,1),...
                indata(tfind,xindim),indata(tfind,xindim)-time_minval,indata(tfind,yindim));
            ksq_time(arr_index,1) = indata(tfind,xindim);
            ksq_time(arr_index,2) = indata(tfind,yindim);
            curr_time = indata(tfind,xindim) + tsample;
            arr_index = arr_index + 1;
        end
    end
end
fclose(fwrite_time);
