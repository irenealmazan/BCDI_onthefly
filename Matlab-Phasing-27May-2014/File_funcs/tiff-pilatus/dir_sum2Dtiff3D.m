name_of_this_file='dir_sum2Dtiff3D';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
this_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

disp('checking directory for .tif files...')


files=rdir([[this_dir],'**/*.tif']);


fnames=char(files.name);

nn=size(fnames);

save_dir=this_dir;

first=fnames(1,:);
last=fnames(end,:);

disp(['first file is - ', first])
disp(['last file is - ',last])
disp('summing....')
sum2Dtiff3D({[first],[last]},1,save_dir);