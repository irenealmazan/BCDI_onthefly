readme for bin_crop_center.m


useage 1:

full_files={'\blah\file1.spe','\blah\file2.spe'}; %full file paths
full_bg={'\blah\bg1.spe','\blah\bg2.spe'}; %full background file paths, OR leave blank, full_bg={};
bin=[1,1] ; % [x,y] Binning.

min_data=300; %threshold for data, everything below this will be removed (i.e set to 0)

aliens=[0]; %use to remove alien scatter, specify [x1,y1,z1,x1,y1,z1] to remove that segment of data, leave as 0 for nothing

nnc=[-20,-20,20,20,10,10]; %does initial padding/cropping. [-20,-20,20,20,10,10] will remove 20 pixels from x (either side), pad 20 to y (either side) and add 10 to z (either side). set nnc=[0] to do nothing. 

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc);

useage 2:

same as above bit now with file_params to pass extra paramters

file_params.subtract_dc = 1;% will subtract the min_data value off the data
rather than just setting everything below this to 0. 

file_params.schot_th = 200; % will do a secondary threshold after binning
useful since some schot noise won't get eliminated by the first threshold
but a second higer one will

file_params.no_center = 1;%, won't center the data, useful for ptycho.
file_params.no_fft_pad=1;%, won't do any fft padding
file_params.pad_ptych=1; %only pad x and y

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc,file_params);
______________________________________________________________________________________________________
SUGGESTED USEAGE:

full_files={'\blah\file1.spe','\blah\file2.spe'}; %full file paths
full_bg={'\blah\bg1.spe','\blah\bg2.spe'}; %full background file paths, OR leave blank, full_bg={};
bin=[1,1] ; % no Binning.

min_data=0; %threshold of 0
aliens=[0]; %no alien removal

nnc=[0]; %no initial padding or cropping 

file_params.subtract_dc = 0;% don't do a DC subtraction

file_params.schot_th = 0; % no secondary thresholding
file_params.no_center = 1;%, don't center the data
file_params.no_fft_pad=0;%, will pad for fft
file_params.pad_ptych=1;% will do x,y only

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc,file_params);
