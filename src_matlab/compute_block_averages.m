function std_arr = compute_block_averages(inparr)

%yields an array of size(length(inparr)/2,3)
%first column is the block size, n
%second is the std (\sigma_n) with a block size, n
%third is BSE = \sigma_n/sqrt(M), where M = length(inparr)/n
%REF: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865156/
format long
len_arr       = length(inparr(:,1));
max_blocksize = floor(len_arr/2); % This number is used to limit the number of blocks. Num_blocks is calcualted wrt total length
std_arr       = zeros(max_blocksize,3); 

fout = fopen('./../../rgtime_data/errcheck.dat','w');

for i = 1:max_blocksize
    block_size = i;
    num_blocks = floor(len_arr/block_size); % to make arrays don't blow
    avg_arr    = zeros(num_blocks,1); % to make sure zeros are not added to array
    cntr = 0;
    fprintf(fout,'blocksize: %d, nsamples: %d\n',i, num_blocks);
    for j = 1:block_size:num_blocks*block_size
        cntr = cntr + 1;
        avg_arr(cntr,1) = mean(inparr(j:j+block_size-1));
        fprintf(fout,'%d\t%d\n',j,avg_arr(cntr,1));
    end
    % brute force calculation - NOT using STD function in MATLAB
    meanval = mean(avg_arr(:,1));
    maxcntr = cntr;
    std_arr(i,1) = block_size;
    net_err = 0;
    for j = 1:maxcntr
        net_err = net_err + (avg_arr(j,1)^2 - meanval^2);
    end
    std_arr(i,2) = sqrt(net_err/(num_blocks-1)); % This is sigma_n (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2865156/)
    std_arr(i,3) = std_arr(i,2)/sqrt(num_blocks);
    clear avg_arr net_err
end