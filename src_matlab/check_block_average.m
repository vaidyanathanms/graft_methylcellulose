%https://stackoverflow.com/questions/12694295/statistical-inefficiency-block-averages
datain = importdata('./../../rgtime_data/check_data.txt');
std_out   = compute_block_averages(datain(:,2));
plot(std_out(:,1),std_out(:,2),'ko')
