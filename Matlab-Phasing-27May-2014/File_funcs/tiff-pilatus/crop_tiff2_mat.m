%converts tif files to matlab (2D or 3D tiff's)
%change the prefix below accrodingly

name_of_this_file='crop_tiff2_mat';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
this_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

disp('checking directory for crop@@@.tif files...')


files=rdir([[this_dir],'crop*.tif']);
fnames=char(files.name);

info = imfinfo(strcat(fnames));
num_images = numel(info);

save_name=fnames(numel(this_dir)+1:end);
save_name=[save_name(1:end-4),'.mat'];

for qq = 1:num_images
    
    ff0 = imread(fnames,qq);
    
    if qq == 0, data=zeros([size(ff0),num_images]);end
    
    data(:,:,qq) = ff0;
    
end

save([this_dir,save_name],'data');
