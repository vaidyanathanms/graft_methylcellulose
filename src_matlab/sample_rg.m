function rgout_data = sample_rg(indata, xindim, yindim, tinarr, tsample, config, gr_mw, sig_val, eps_val)

format long
tol = 1e-7;
%% Sample uniform data from input
inp_tminval  = min(tinarr); inp_tmaxval = max(tinarr); 
out_nsamples = (inp_tmaxval-inp_tminval)/tsample;
rgout_data   = zeros(floor(out_nsamples),2);

% store sampled values in a separate file
fylename = sprintf('../../rgtime_data/sampled_data/rgsampletime_conf_%d_grmw_%d_sig_%g_eps_%g',...
    config,gr_mw,sig_val,eps_val);
fwrite_time = fopen(fylename,'w');
fprintf(fwrite_time,'%s\t %s\t %s\t %s %g\n','Step (Unshifted)','Unshifted time','Shifted time','Kappasq with sample time: ',tsample);

% Find minimum value and then start from that point while subtracting the
% min time value (DO NOT subtract the step value)
mintime_flag = -1;
for tfind = 1:length(indata(:,xindim))
    if indata(tfind,xindim) >= inp_tminval
        tmin_index = tfind;
        time_minval = indata(tfind,xindim);
        fprintf(fwrite_time,'%g\t%g\t%g\t%g\n',indata(tfind,1),...
            indata(tfind,xindim),indata(tfind,xindim)-time_minval,indata(tfind,yindim));
        mintime_flag = 1;
        break;
    end
end

if mintime_flag == -1
    fprintf('Simulation has not reached the minimum time prescribed \n');
    fprintf('Current time: %g/Min time %g \n',indata(tfind-1,xindim),inp_tminval)
    return;
end

curr_time       = time_minval;
rgout_data(1,1) = curr_time;
rgout_data(1,2) = indata(tfind,yindim);
arr_index       = 2;
curr_time       = indata(tfind,xindim) + tsample;

for tfind = tmin_index+1:length(indata(:,xindim))
    if indata(tfind,xindim) <= inp_tmaxval %% Sample for values less than max value
        if abs(curr_time-indata(tfind,xindim)) < tol %% Sample only at equispaced intervals
            fprintf(fwrite_time,'%g\t%g\t%g\t%g\n',indata(tfind,1),...
                indata(tfind,xindim),indata(tfind,xindim)-time_minval,indata(tfind,yindim));
            rgout_data(arr_index,1) = indata(tfind,xindim);
            rgout_data(arr_index,2) = indata(tfind,yindim);
            curr_time = indata(tfind,xindim) + tsample;
            arr_index = arr_index + 1;
        end
    end
end

fclose(fwrite_time);
