function process_and_write(directory_name,fylelist,fylecolumns,outfyle)
fout = fopen(outfyle,'w');
cd(directory_name)
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
    for lskip = 1:fylecolumns(1,4)
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end
        if i == 1
            spl_tline = strsplit(tline);
            fprintf(fout,'%s\t%s\n',spl_tline{fylecolumns(1,2)},spl_tline{fylecolumns(1,3)});
        end
    end    
    while true
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end    
        spl_tline = strsplit(tline);
        xval = str2double(spl_tline{fylecolumns(1,2)});
        yval = str2double(spl_tline{fylecolumns(1,3)});
        fprintf(fout,'%g\t%g\n',xval,yval);
    end
    fclose(finp);
end
fclose(fout);
        
