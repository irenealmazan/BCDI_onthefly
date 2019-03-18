%the base dir where the files are
data_dir='/Users/jesseclark/Documents/MATLAB/data_analysis/APS-Ptych-July2012/';            

%the dir to save the processed files
base_dir='/Users/jesseclark/Documents/MATLAB/data_analysis/APS-Ptych-July2012/Processed/';

%scan numbers to process
scan_numbers=[75,77];

%the size of the number as a string. 075 is 3, 75 is 2 and 7 is 1
num_str_size=3;

%file prefix's
prefix='uclpup12-2b_*';

%prefix for the scan
scan_prefix='S';

%file type
ftype='.tif';

for qq=1:numel(scan_numbers)
        
    num=check_numlength(scan_numbers(qq),num_str_size);

    %dir where the files to be added together are
    cur_dir=[data_dir,scan_prefix,num,'/'];
        
    %create save directory
    save_dir=[base_dir,scan_prefix,num,'/'];%,

    %create the directory if not already made
    if isdir([save_dir]) ~= 1,mkdir([save_dir]);end

    fnames=rdir([cur_dir,prefix,ftype]);
    
    first=fnames(1).name;
    last=fnames(end).name;

    n1=extract_number_from_string(first);
    n2=extract_number_from_string(last);
    
    matname=[scan_prefix,num,'-',num2str(n1),'-',num2str(n2),'.mat'];
 
    %use this one if the files should be created within this func
    %data=Tiff2D_2_Tiff3D_orientV2({first,last},1,[save_dir,scan_prefix,num,'-',],'APS-pil');
    
    %use this one if you want to just pass through the files returned by
    %rdir
    data=Tiff2D_2_Tiff3D_orientV2(fnames,1,[save_dir,scan_prefix,num,'-',]);
    
    
    save([save_dir,matname],'data')

end
    

1;

 %copy matlab phasing with correct data dir etc
% template_dir='/Volumes/HENRIQUE/Users/jesseclark/Documents/MATLAB/data_analysis/DLS-I07-Au-July2012/ZTemplate-DO-NOT-REMOVE/';
% string=['cp ',template_dir,'Matlabphasing_ver1_1.m ',save_dir,'/rec/'];
% system(string)
% %now change the data dir and data file to be the correct one
% mphasefile=[save_dir,'/rec/Matlabphasing_ver1_1.m'];
% modify_file(mphasefile,'params.data_dir=''''',['params.data_dir=''',save_dir,'/data/','''']);
% modify_file(mphasefile,'params.files={''''}',['params.files={''',matname,'''}']);
