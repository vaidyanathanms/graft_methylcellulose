function outdata = compute_autocorr(datain, in_dim)
%compute spectral product
fprintf('Computing Autocorrelation Function ..\n');
nframes = length(datain);
outdata = zeros(nframes,1);

for tinc = 0:nframes-1
    ifin = nframes - tinc;
    
    for i = 1:ifin
   
        tim = i + tinc;
        outdata(tinc+1) = outdata(tinc+1) + (datain(tim,in_dim)*datain(i,in_dim));
        
    end
    
    outdata(tinc+1) = outdata(tinc+1)/(ifin);
    
end
