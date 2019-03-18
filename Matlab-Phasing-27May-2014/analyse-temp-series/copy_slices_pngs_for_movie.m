function copy_slices_pngs_for_movie

amp_ph='Ph';

orient='xz'; %xy,xz,zy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_of_this_file='copy_slices_pngs_for_movie';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
fdir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

fstring=['*I-',amp_ph,'-',orient,'.png'];

fullnames=rdir([fdir,'**/*',fstring]);

nnames=size(fullnames,1);

save_dir=[fdir,'movies/',amp_ph,'-',orient,'/'];

if isdir(save_dir) ~= 1,mkdir(save_dir);end

for qq=1:nnames
    
    string=['cp ',char(fullnames(qq).name),' ',save_dir];
    system(string)
    
end



end