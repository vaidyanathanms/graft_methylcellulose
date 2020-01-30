function process_and_write_distCOM(directory_name,fylelist,fylecolumns,outfyle,nchains)
fout = fopen(outfyle,'w');
cd(directory_name)
%open all files within the fylelist
nfyles = numel(fylelist); %number of files of the type

if nchains ~=2
    error('Functionality for #chains not equal to 2 is not implemented');
end

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
            fprintf(fout,'%s\t%s\t%s\t%s\t%s\n','Timestep','dx','dy','dz','dCOM');
        end
    end    
    while true
        tline = fgetl(finp);
        if ~ischar(tline)
            break;
        end    
        spl_tline = strsplit(tline);
        if size(spl_tline) ~= 5
            fprintf('Not all details are present in %\n',tline)
            error('Exiting at file: %s\n',fylename)
        end
        
        tval  = str2double(spl_tline{1});
        xdist = str2double(spl_tline{2});
        ydist = str2double(spl_tline{3});
        zdist = str2double(spl_tline{4});
        dcom  = str2double(spl_tline{5});
        fprintf(fout,'%g\t%g\t%g\t%g\t%g\n',tval,xdist,ydist,zdist,dcom);
    end
    fclose(finp);
end
fclose(fout);
        