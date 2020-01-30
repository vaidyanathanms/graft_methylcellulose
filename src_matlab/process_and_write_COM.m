function process_and_write_COM(directory_name,fylelist,fylecolumns,outfyle,nchains)
fout = fopen(outfyle,'w');
cd(directory_name)
%open all files within the fylelist
nfyles = numel(fylelist); %number of files of the type
com_vals = zeros(nchains,3); box_vals = zeros(3);


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
        if size(spl_tline) ~= 8
            fprintf('Not all details are present in %\n',tline)
            error('Exiting at file: %s\n',fylename)
        end
        
        chid = str2double(spl_tline{2});
        com_vals(chid,1) = str2double(spl_tline{3});
        com_vals(chid,2) = str2double(spl_tline{4});
        com_vals(chid,3) = str2double(spl_tline{5});
        
        if chid == 1
            tval = str2double(spl_tline{1});
            box_vals(1) = str2double(spl_tline{6});
            box_vals(2) = str2double(spl_tline{7});
            box_vals(3) = str2double(spl_tline{8});
        end
        
        if chid == nchains
            xdist = com_vals(chid,1)-com_vals(chid-1,1);
            ydist = com_vals(chid,2)-com_vals(chid-1,2);
            zdist = com_vals(chid,3)-com_vals(chid-1,3);
            xdist = xdist - box_vals(1)*round(xdist/box_vals(1));
            ydist = ydist - box_vals(2)*round(ydist/box_vals(2));
            zdist = zdist - box_vals(3)*round(zdist/box_vals(3));
            dcom  = sqrt(xdist^2+ydist^2+zdist^2);
            fprintf(fout,'%g\t%g\t%g\t%g\t%g\n',tval,xdist,ydist,zdist,dcom);
            com_vals = zeros(nchains,3); box_vals = zeros(3);
        end
    end
    fclose(finp);
end
fclose(fout);
        
