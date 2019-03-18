function data2vtk(dir_file)
%jclark
%load data from either a .rec file or params file and 
%output as vtk


%load the params file
pfile=rdir([dir_file,'**/*PARAMS.mat']);    
load(pfile(1).name)

%
data_dir0=params.data_dir;
nnc=params.nnc;
files=params.files;
back=params.back;
min_data=params.min_data;
params.no_hist=1;
aliens=params.aliens;
bin=params.binning;

full_files=strcat(data_dir0,files);
if numel(back) ~= 0,full_bg=strcat(data_dir0,back);else full_bg=[];end

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc,params);

filename=[dir_file,'data.vtk'];
savevtk2scalar(data,filename,[],1)

end

