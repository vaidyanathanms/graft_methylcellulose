function [min_indexcut,max_indexcut] = find_general_cutoff_minmax(data,mincut,maxcut)

lenval = length(data(:,1));
min_indexcut = floor(0.8*lenval);
max_indexcut = lenval;

%find minimum cut
for i = 1:lenval
    if data(i,1) > mincut
        min_indexcut = i;
        break;
    end
end

%find maximum cut
for i = 1:lenval
    if data(i,1) > maxcut
        max_indexcut = i;
        break;
    end
end