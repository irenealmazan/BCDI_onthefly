%using Andrejs stuff to figure out where the highest counts are


load('/Users/a4894z/Dropbox/Battery project/RockingData.mat')

plot(data.run_arr,data.P(1,:),'r')

%indices don't quite line up so have to find the new ones
beg_num=868;

for i=1:20
	ind(i) = find(data.run_arr==beg_num+i-1);
	beg_num+i-1
end

%during phase transition
snums = 886:896;

for i=1:length(snums)
	ind = find(data.run_arr==snums(i));
	disp('for scan number')
	snums(i)
	disp('total counts')
	data.P(1,ind)
	
end