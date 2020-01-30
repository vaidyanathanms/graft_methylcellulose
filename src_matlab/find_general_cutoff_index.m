function indexcut = find_general_cutoff_index(data,cutoffval)

lenval = length(data(:,1));
indexcut = floor(0.8*length(data(:,1)));
for i = 1:lenval
    if data(i,1) > cutoffval
        indexcut = i;
        break;
    end
end