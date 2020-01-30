function indexcut = find_cutindex(data,cutoffval)

lenval = length(data(:,1));
indexcut = length(data(:,1));
for i = 2:lenval
    if data(i,1) > data(i-1,1) && data(i,1) < 0.3 % initially the correlation will increase. So should account for that
        indexcut = i-2;
    elseif data(i,1) < cutoffval
        indexcut = floor(0.9*i);
        break;
    end
end