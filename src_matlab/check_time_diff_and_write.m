function check_time_diff_and_write(directory_name,fylelist,fylecolumns,outfyle,tdiff)

%open file pointers to write all chain details separately to plot
fout = fopen(outfyle,'w');

cd(directory_name);
%open all files within the fylelist
nfyles = numel(fylelist); %number of files of the type

for i = 1:nfyles
    fylename = fylelist(i).name;
    if exist(fylename,'file') ~= 2
        fprintf('%s does not exist\n',fylename)
        continue;
    end
        
    all_data = importdata(fylename);
    tdata    = all_data.data(:,1);
    
    %proceed if and only if the time difference matches the input time
    %diff.
    if tdata(2,1) - tdata(1,1) ~= tdiff
        continue;
    end
  
    fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Timestep','ChainID',...
        'max_eig','Rg^2','asphericity','acylindricity','shapefactor',...
        'prolateness');
   
    
    finp = fopen(fylename,'r');
    fprintf('File under process %s\n',fylename);
    for lskip = 1:fylecolumns(1,2) %headerlines
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end
        
    end
    
    while true
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end    
        spl_tline = strsplit(tline);
        if size(spl_tline) ~= 5
            fprintf('Not all details are present in %s \n',tline)
            error('Exiting at file: %s\n',fylename)
        end

        tval = str2double(spl_tline{1});
        chid = str2double(spl_tline{2});
        lam_xsq = str2double(spl_tline{3});
        lam_ysq = str2double(spl_tline{4});
        lam_zsq = str2double(spl_tline{5});
        sum_lam = lam_xsq + lam_ysq + lam_zsq; %rgsq
        asphere = lam_zsq - 0.5*(lam_xsq + lam_ysq);
        acylndr = lam_ysq - lam_xsq;
        shapefac = 1.5*((lam_xsq^2 + lam_ysq^2 + lam_zsq^2)/(sum_lam^2)) - 0.5;
        norm_asphere = asphere/sum_lam;
        
        %compute prolateness
        
        mean_lam = sum_lam/3;
        prod_lam = (lam_xsq - mean_lam)*(lam_ysq - mean_lam)*(lam_zsq - mean_lam);
        num_val  = 4*prod_lam;
        lam_dev  = (lam_xsq - mean_lam)^2 + (lam_ysq - mean_lam)^2 + (lam_zsq - mean_lam)^2;
        den_val  = ((2/3)*lam_dev)^(1.5);
        prolate  = num_val/den_val; 
        
        fprintf(fout,'%d\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',tval,chid,...
            lam_zsq,sum_lam,asphere,acylndr,shapefac,prolate);
    end
    fclose(finp);
end
fclose(fout);



