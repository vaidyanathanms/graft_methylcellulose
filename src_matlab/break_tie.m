function ident = break_tie(sim_data,shapedata)
%% Rules: if equal probability is found, return the identity of each cases concatenated

% check if NA is present as the major weight. If so exit with error
indices = (find(sim_data(:,1) == max(sim_data(:,1))));
for cntr = 1:length(indices)
    if strcmp(shapedata{indices(cntr)}{1},'NA')
        fprintf('One of the tied configurations is "NA": Run more simulations\n')
    end
end


% possible combinations: TDT, DTSH, DTML
% return the identity of the string for equal matches
if length(indices) ~=2
    fprintf('Error: Can work with only two equal probability cases: %d\n',length(indices));
    return;
end
ind_sort = sort(indices); % to make sure that the order of combinations fall in the possible combinations given above
ident = strcat(shapedata{ind_sort(1)}{1},shapedata{ind_sort(2)}{1});



