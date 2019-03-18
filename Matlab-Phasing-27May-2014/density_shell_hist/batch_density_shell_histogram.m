function batch_density_shell_histogram
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

name_of_this_file='batch_density_shell_histogram';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
script='CVL-AMP.rec';
n_sc=numel(script);
files=rdir([[dir],'**/*',script]);
fnames=char(files.name);

nn=size(fnames,1);

support=[];
shells=5;

for qq=3:nn
    
    save_dir=strtrim(fnames(qq,:));
    save_dir=save_dir(1:end-numel(script));
    load(strtrim(fnames(qq,:)),'-mat')
    density_shells_histogram(array,support,shells,save_dir )
    
end

end

