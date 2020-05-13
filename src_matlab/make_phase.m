%% To generate phase diagram from the snapshots
% Data are sorted in final_configs.xlsx
% Sheet!all_data in final_configs.xlsx contains consolidated data

clc;
clear;
close all;
format long;

%% Convertors - Add more structures in shape_define

% Color Codes for Plot
green = [0 0.5 0.0]; gold = [0.9 0.75 0]; orange = [0.91 0.41 0.17]; brown=[0.6 0.2 0];

% Toroid - 'T' = 1
% Distorted Toroid/Compact Hairpins/Eight-likw Toroids - 'DT'/'CH'/'ET' = 2
% Multiloop/Sphere-like - 'ML'/'SL' = 4
% Stretched Hairpins - 'SH' = 6
% Open Structures - 'O' = 7
% Random Intermediates - 'RI' = 0
% Not available - 'NA' = -1
% Not equilibrated - 'NEQ' = -1

shape_define = {{'T',1,'k','d'}; %{name,shape_id,color,markertype}
    {'DT',2,orange,'^'};
    {'CH',2,orange,'^'};
    {'ET',2,orange,'^'};
    {'ML',3,'m','h'};
    {'SL',3,'m','h'};
    {'SH',4,gold,'s'};
    {'RI',4,gold,'s'};
    {'O',5,green,'o'};
    {'NA',-1,'none','none'}
    {'NEQ',-1,'none','none'}};
num_shapes = length(shape_define);

%equal_prob_cases - vertices for plotting (see break_tie.m and plotCustMarkDemo.m)
eq_prob = {{'TDT',6,[0 -0.25 0],[0.5 0 -0.5],[0 0.33 0],[0.5 -0.5 -0.5],'k',orange};
    {'DTSH',7,[0 -0.25 0],[0.5 -0.5 -0.5],[0 0.25 0.25 0],[0.5 0.5 -0.5 -0.5],orange,gold};
    {'DTML',8,[0 -0.25 0],[0.5 0 -0.5],[0 0.125 0.25 0.125 0.25 0.125 0],[0.5 0.25 0.25 0 -0.25 -0.25 -0.5],orange,'m'}}; %{name,id,vertx1,verty1,vertx2,verty2,col1,col2}

%% Change this if new shapes are added
legendchar{1} = 'T, Toroid';
legendchar{2} = 'DT, Distorted Toroid';
legendchar{3} = 'SL, Sphere-like/Multi-loop';
legendchar{4} = 'SH, Stretched Hairpins';
legendchar{5} = 'OR, Open-loop/Rod-like';
% legendchar{6} = 'T/DT, Toroid/Distorted Toroid';
% legendchar{7} = 'DT/SH, Distorted Toroid/Stretched Hairpins';
% legendchar{8} = 'DT/ML, Distorted Toroid/Multi-loop';
max_flag_legends = length(legendchar);

%% Input data

eps_arr = [0.8,1,1.2];
sig_arr = [0,0.01,0.05,0.1,0.15,0.2,0.25,0.3];

shape_arr = cell(length(eps_arr),length(sig_arr));  %name of the shape
prob_arr  = zeros(length(eps_arr),length(sig_arr)); %probability of finding that shape
flag_idequal_prob = zeros(length(eps_arr),length(sig_arr)); %flag for equal probability cases

fyle_id = fopen('../../snapshots/all_final_configs.dat','r');
if fyle_id <= 0
    error('File does not exist \n');
end

outfyle_id = fopen('../../snapshots/most_probable_config.txt','w');
fprintf(outfyle_id,'%s\t%s\t%s\t%s\t%s\n','epsilon','sigma','config','maxruns','probability');

find_eps = -1; sig_header = -1;

while ~feof(fyle_id)
    
    tline = fgetl(fyle_id); %get line
    
    if isempty(strtrim(tline))
        continue; %skip empty lines
    end
    
    all_strings = strsplit(strtrim(tline)); %split line
    
    len_strings = length(all_strings); %find number of words in the line
    
    if strcmp(all_strings{1},'eps') % eps line
        
        epsval = str2double(all_strings{2});
        
        %check whether it is present in the main eps_array
        find_eps = -1;
        for eps_cnt = 1:length(eps_arr)
            if epsval == eps_arr(eps_cnt)
                find_eps = 1;
                eps_num  = eps_cnt; %index corresponding to the eps_arr
                break;
            end
        end
        if find_eps == -1
            error('Unknown epsilon value in the file: %g',eps_val);
        end
        sig_header = -1; %check for sigma header for this epsilon value
        fprintf('Analyzing eps = %g \n', epsval);
        
    elseif strcmp(all_strings{1},'sigma/conf') %sigma/conf line
        
        if find_eps ~= 1
            error('No parent epsilon value found in the file');
        end
        
        sig_header = 1; %found sigma header
        
        
    else %process main lines and loop them
        
        if sig_header ~= 1 || find_eps ~= 1
            error('No sigma/epsilon header found');
        end
        
        sigval = str2double(all_strings{1});
        
        %check whether it is present in the main array
        find_sig = -1;
        for sig_cnt = 1:length(sig_arr)
            if sigval == sig_arr(sig_cnt)
                find_sig = 1;
                sig_num  = sig_cnt; %index corresponding to the sig_arr
                break;
            end
        end
        if find_sig == -1
            error('Unknown sigma value in the file: %g',sigval);
        end
        
        num_configs = len_strings-1; %total configuration per sigma
        score_arr_shapes = zeros(num_shapes,1); %zero all counter arrays for this eps/sig value
        num_NA_vals      = 0; %count the number of NA vals
        
        for conf_cnt = 1:num_configs %now for each sigma loop through all the columns (num_configs)
            config_string = all_strings{conf_cnt+1}; % +1 counter because first element is sigma value
            %check whether it is present in the set of configurations at the beginning
            find_shp = -1;
            for shp_cnt = 1:length(shape_define)
                if strcmp(config_string,shape_define{shp_cnt}{1})
                    find_shp = 1;
                    shp_num  = shp_cnt; %index of the shape in the shape_define array
                    if strcmp(config_string,'NA')
                        num_NA_vals = num_NA_vals + 1; %count this separately for finding the net probability of a given shape
                    end
                    break;
                end
            end
            if find_shp == -1
                error('Unknown sigma value in the file %s\n',config_string);
            end
            
            %Update the counter for the particular configuration
            score_arr_shapes(shp_num,1) = score_arr_shapes(shp_num,1) + 1;
        end
        
        
        %% Find most probable configuration for the given sigma/eps
        %new edit: add probability of find the configuration: apr-30-2020
        
        %If there exists tie between configurations' maximum score
        if length(find(score_arr_shapes(:,1) == max(score_arr_shapes(:,1)))) == 1
            [ncntsval,index] = max(score_arr_shapes); %format [m,i] = max(A)
            shape_arr{eps_num,sig_num} = shape_define{index}{1};
        else %break tie
            fprintf('Tie for  eps/sig: %g\t%g\n', epsval,sigval);
            shape_arr{eps_num,sig_num} = break_tie(score_arr_shapes, shape_define); %see rules in break_tie
            flag_idequal_prob(eps_num,sig_num) = 1;
            indices_out = find(score_arr_shapes(:,1) == max(score_arr_shapes(:,1)));
            index = min(indices_out); %this is to get the index of the "minimum" of the two configurations
            [ncntsval,minindex] = max(score_arr_shapes(:,1)); %this is to count the number of repetitions
        end
        
        %if NA is the best - write output as error - need more simulations
        if strcmp(shape_arr{eps_num,sig_num},'NA')
            fprintf('Cannot plot (NA is major config) for eps/sig: %g\t%g\n', epsval,sigval);
        end
        
        prob_arr(eps_num,sig_num)  = ncntsval/(num_configs-num_NA_vals); %probability of finding that configuration
        fprintf(outfyle_id,'%g\t %g\t %s\t %d\t %g\n',epsval,sigval,shape_arr{eps_num,sig_num},...
            num_configs-num_NA_vals,prob_arr(eps_num,sig_num));
        
    end %end all if-cases
    
end % end tline loop


%% Plot all data

hphase = figure;
hold on
box on
set(gca,'FontSize',16)
xlabel('$\Sigma$','FontSize',24,'Interpreter','Latex')
ylabel('$\epsilon_{gg}$','FontSize',24,'Interpreter','Latex')

point_pref = 14;
flag_legend = zeros(max_flag_legends,1);
ax_cntr = 0; flagshp_cntr = 0;

for i = 1:length(eps_arr)
    for j = 1:length(sig_arr)
        if isempty(cell2mat(shape_arr(i,j)))
            continue;
        end
        config_string = shape_arr(i,j);
        find_shp = -1;
        if flag_idequal_prob(i,j) == 0 %for NOT equal probability cases
            for shp_cnt = 1:length(shape_define)
                if strcmp(config_string,shape_define{shp_cnt}{1})
                    find_shp = 1;
                    shp_num  = shp_cnt; %index of the shape in the shape_define array
                    break;
                end
            end
            if find_shp == -1
                error('Unknown shape value while plotting %s\n',cell2mat(shape_arr(i,j)));
            end
            
            ax_cntr = ax_cntr + 1;
            ax(ax_cntr) = plot(sig_arr(j),eps_arr(i),'Color',shape_define{shp_num}{3},'Marker',...
                shape_define{shp_num}{4},'MarkerFaceColor',shape_define{shp_num}{3},...
                'MarkerSize',point_pref*prob_arr(i,j),'LineStyle','None');
            flag_cntr = shape_define{shp_cnt}{2};
        
            if flag_legend(flag_cntr,1) == 0
                flag_legend(flag_cntr,1) = ax_cntr;
            end
        else %for equal probability cases
            for shp_cnt = 1:length(eq_prob)
                if strcmp(config_string,eq_prob{shp_cnt}{1})
                    find_shp = 1;
                    shp_num  = shp_cnt; %index of the shape in the equal_prob_cases array
                    break;
                end
            end
            if find_shp == -1
                error('Unknown shape value while plotting %s\n',cell2mat(shape_arr(i,j)));
            end
            
            % myCustMarker(x,y,vert_thex_X1,vert_thex_Y1,vert_thex_X2,vert_thex_Y2,scaleval,mf1,me1,mf2,me2)
            
            ax_cntr = ax_cntr + 1;
            ax(ax_cntr) = myCustMarker(sig_arr(j),eps_arr(i),...
                eq_prob{shp_num}{3},eq_prob{shp_num}{4},eq_prob{shp_num}{5},eq_prob{shp_num}{6},...
                0.005*point_pref*prob_arr(i,j),eq_prob{shp_num}{7},eq_prob{shp_num}{7},eq_prob{shp_num}{8},eq_prob{shp_num}{8},hphase);
            flag_cntr = eq_prob{shp_cnt}{2};
            
            % To represent legend only once
        end
        
        
        
    end
end

%Represent legend only once
legendarr   = zeros(max_flag_legends,1);
for i = 1:max_flag_legends
    legendarr(i,1)   = ax(flag_legend(i));
end

legend(legendarr,legendchar,'FontSize',12,'Location','Northeast')

ylim([0.75 1.52]);
yticks([0.8 1.0 1.2]);
saveas(hphase,'./../../all_figures/fig_phasedia.png');



