function align_SPE2tif(dir_file,fft_pad)

%converts a series of spe's to a tif.  makes it easier for alien removal
%in someething like imagej.  assumes that matlabphaing has been run.  will
%use the params file for various inputs (including file names).
%
% i think the preferred usage is align_SPE2tif([pwd,'/'],1)
% this will do an intial padding then crop to the original size
% (params.return_orig_size=1).  This means that there is minimal spill over
% if the pattern is not so centered

%%
try
    data_dir=dir_file;
catch
    name_of_this_file='align_SPE2tif';
    dir_file=which(name_of_this_file);
    dir_file=dir_file(1:findstr(dir_file,name_of_this_file)-1);
    data_dir=dir_file;
end


%pfile=rdir([dir_file,'rand-starts/','**/*PARAMS.mat']);
%if numel(pfile) == 0
%    pfile=rdir([dir_file,'**/*PARAMS.mat']);
%else

pfile=rdir([dir_file,'**/*PARAMS.mat']);    

counter=0;
qq=0;

while counter == 0      %searches for params files that used the .spe's
    qq=1+qq;
    temp=pfile(qq).name;
    load(temp);
    temp_n=char(params.files(1));
    counter=strcmp(lower(temp_n(end-3:end)),'.spe');
    if qq == max(size(pfile)),counter = 1;end           %incase it doesn't find one
end
 
load(temp)

%pfile=pfile(1).name;
%load(pfile)

%bin=params.binning; %get the binning, want to display with minimal binning though
bin=[1,1];

% try 
%     aliens=params.aliens;
% catch
%     aliens=[];
% end
aliens=[0];

files=params.files;
back=params.back;
min_data=params.min_data;
params.no_hist=1;

try
    fft_pad;
catch
    fft_pad=1;
end

try
    data_dir0=params.data_dir;
catch
    data_dir0=data_dir;
end

if fft_pad == 1,params.no_fft_pad=0;else params.no_fft_pad=1;end

params.schot_th=0;

try
    nnc=params.nnc;             % not support yet, will be used for initial cropping
catch
    nnc=[0];
end
nnc=0;

params.return_orig_size=1;
params.do_2D=0;

full_files=strcat(data_dir0,files);
if numel(back) ~= 0,full_bg=strcat(data_dir0,back);else full_bg=[];end

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc,params);

%get the number for the file name
nfiles=numel(files);
sname=[];
for qq=1:nfiles

    a=cell2mat(files(qq));
    numb=num2str(sscanf(a(strfind(a,'-')+1:numel(a)),'%i'));

    if numel(numb) == 0,numb=num2str(sscanf(a(strfind(a,'_')+1:numel(a)),'%i'));end

    numbs(qq)=str2num(numb);
    
    szn=size(numb);
    if szn(1) > 1,numb=numb(end,:);end
    
    if qq ~= nfiles, sname=[sname,numb,'-'];else sname=[sname,numb];end

    
end
    
mat2tif(data,[dir_file,'/Combined-',sname,'.tif'])

end



