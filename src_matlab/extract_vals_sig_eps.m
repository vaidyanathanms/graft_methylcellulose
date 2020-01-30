function out_arr = extract_vals_sig_eps(main_data,sigcol,epscol,datacol,eps_dumarr,sigma_dumarr)
out_arr = zeros(length(eps_dumarr),length(sigma_dumarr));
sig_allvals  = main_data.data(:,sigcol);
eps_allvals  = main_data.data(:,epscol);
main_allvals = main_data.data(:,datacol);
len_fyle = length(main_allvals(:,1));

for fcnt = 1:len_fyle
    sigval = sig_allvals(fcnt,1);
    epsval = eps_allvals(fcnt,1);
    
    epsindex_flag = -1;
    for epcnt = 1:length(eps_dumarr)
        if eps_dumarr(epcnt) == epsval
            epsindex = epcnt;
            epsindex_flag = 1;
            break;
        end
    end
    
    if epsindex_flag == -1
        fprintf('Did not find epsilon value = %s at grMW = %d\n',...
            eps_dumarr(epcnt), gr_mw);
        continue;
    end
    
    sigindex_flag = -1;
    
    for sigcnt = 1:length(sigma_dumarr)
        if sigma_dumarr(sigcnt) == sigval
            sigindex = sigcnt;
            sigindex_flag = 1;
            break;
        end
    end
    
    if sigindex_flag == -1
        fprintf('Did not find sigma value = %s at grMW = %d\n', ...
            sigma_dumarr(sigcnt), gr_mw);
        continue;
    end
    
    out_arr(epsindex,sigindex) = main_allvals(fcnt);
    
end