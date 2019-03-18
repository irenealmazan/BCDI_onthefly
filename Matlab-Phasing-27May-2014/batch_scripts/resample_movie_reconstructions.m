function  resample_movie_reconstructions(dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
scale=0.75;

try
    dir;
catch
    name_of_this_file='resample_movie_reconstructions';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
end

script='.rec';

files=rdir([[dir],'**/*',script]);
fnames=char(files.name);

for qq=1:size(fnames,1)
    
    name_qq=strtrim(fnames(qq,:));
    disp(name_qq(length(dir)+1:end))
    load(name_qq,'-mat');
    Farr=imresize(ifftshift(fftn(fftshift(array))),scale,'nearest');
    array=real(ifftshift(ifftn(fftshift(Farr))));
    save(name_qq,'array');
    
    
end



end

