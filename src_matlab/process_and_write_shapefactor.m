function process_and_write_shapefactor(directory_name,fylelist,fylecolumns,outfyle,numchains)

%open file pointers to write all chain details separately to plot
fout = fopen(outfyle,'w');
fptr_out = zeros(numchains,1);
for i = 1:numchains
    fyle_split = strsplit(outfyle,'.dat');
    fchfyle = strcat(fyle_split{1},'_',num2str(i),'.dat');
    fptr_out(i,1) = fopen(fchfyle,'w');
    fprintf(fptr_out(i,1),'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Timestep','ChainID',...
        'max_eig','Rg^2','asphericity','acylindricity','shapefactor','prolateness');
end

cd(directory_name);
%open all files within the fylelist
nfyles = numel(fylelist); %number of files of the type

for i = 1:nfyles
    fylename = fylelist(i).name;
    if exist(fylename,'file') ~= 2
        fprintf('%s does not exist\n',fylename)
        continue;
    end
    finp = fopen(fylename,'r');
    fprintf('File under process %s\n',fylename);
    for lskip = 1:fylecolumns(1,2) %headerlines
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end
        if i == 1
            fprintf(fout,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Timestep','ChainID',...
                'max_eig','Rg^2','asphericity','acylindricity','shapefactor',...
                'prolateness');
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
        
        fprintf(fout,'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',tval,chid,...
            lam_zsq,sum_lam,asphere,acylndr,shapefac,prolate);
        fprintf(fptr_out(chid,1),'%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',...
            tval,chid,lam_zsq,sum_lam,asphere,acylndr,shapefac,prolate);
    end
    fclose(finp);
end
fclose(fout);

for i = 1:numchains
    fclose(fptr_out(i,1));
end


