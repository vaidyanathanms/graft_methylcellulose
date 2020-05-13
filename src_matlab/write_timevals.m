function write_timevals(simdirname,job_fylelist,outfname,outf2name)
%% Extract time out of job files

% Open output file and change directory to working directory
foutid = fopen(outfname,'w'); %detailed time-dt values
fsumid = fopen(outf2name,'w'); %summary of dt values
fprintf(foutid,'%s\t %s\n','Timestep','delta t');
fprintf(fsumid,'%s\t%s\t%s\n','Tstep-begin','Tstep-end','delta t');
cd(simdirname);
dtval = 0.005; %default value -- set default at beginning - after that it will be the previous value

for i = 1:length(job_fylelist)
    if contains(job_fylelist(i).name,'ana') %skip all job*ana* files
        fprintf('Skipping analysis file %s\n', job_fylelist(i).name)
        continue
    end
    
    fin_name = dir(job_fylelist(i).name); %check if file is empty
    if fin_name.bytes == 0
        fprintf('No data in %s\n',job_fylelist(i).name);
        continue
    end
    
    fprintf('Analyzing %s\n',job_fylelist(i).name);
    freadid = fopen(job_fylelist(i).name,'r');
    
    while ~feof(freadid)
        strarr = strsplit(strtrim(fgetl(freadid))); 
        if strcmp('timestep',strarr{1}) % check if there is timestep command
            dtval = str2double(strarr{2});
        elseif strcmp('Step',strarr{1}) % check outputs where Steps are output
            strarr = strsplit(strtrim(fgetl(freadid))); sumwrite = -1;
            while ~feof(freadid) && isnumeric(str2double(strarr{1})) && ~isnan(str2double(strarr{1}))
                fprintf(foutid,'%g\t%g\n',str2double(strarr{1}),dtval); %write everything to detailed time-dt file
                if sumwrite == -1
                    fprintf(fsumid,'%g\t',str2double(strarr{1})); %write the beginning timestep to summary file
                    sumwrite = 1;
                end
                lasttimestep = str2double(strarr{1});
                strarr = strsplit(strtrim(fgetl(freadid)));
            end
            fprintf(fsumid,'%g\t%g\n',lasttimestep,dtval);
        end
    end
    fclose(freadid);
end
fclose(foutid);
fclose(fsumid);

