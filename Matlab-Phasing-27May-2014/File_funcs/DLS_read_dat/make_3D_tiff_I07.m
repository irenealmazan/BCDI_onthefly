function [ output_args ] = make_3D_tiff_I07(dat_numb,data_dir,save_dir)
%jclark
if exist('data_dir') ~= 1
    data_dir='/Volumes/HENRIQUE/Users/jesseclark/Data/DLS-I07-AuTe-Oct2013/';
end

if exist('save_dir') ~= 1
    save_dir='/Volumes/HENRIQUE/Users/jesseclark/Data/DLS-I07-AuTe-Oct2013/Processed/';
end

if isdir(extract_dir_from_string(save_dir)) ~= 1,mkdir(extract_dir_from_string(save_dir));end

dpref='pilatus3/';
fpref='p3Image';

[first last im_dir]=DLS_extract_fnums([data_dir,dat_numb],fpref);

if isempty(im_dir) ~=1,dpref=[im_dir,'/'];end

first=[data_dir,dpref,first];
last=[data_dir,dpref,last];

[ data ] = Tiff2D_2_Tiff3D_orientV2({first,last},1,[save_dir]);



end


function [ data ] = Tiff2D_2_Tiff3D_orientV2( files,accum,save_dir,orient)
%jclark
%function Tiff2D_2_Tiff3D( files,accum,save_dir) old usage jnc july2012
%takes a series of 2D tiff's and saves as a 3d tiff.
%use sum2Dtiff3D.m to save as .mat file.
%specify the first and last files eg. files={first_file,last_file}
%specify more than 2 files to load them by filename
%use accum to specify the accumulations

if exist('orient') ~= 1,orient=[];end
    

if isstruct(files) %for use when rdir has been used
    
    first = num2str(extract_number_from_string(files(1).name));
    last = num2str(extract_number_from_string(files(end).name));
    
    name=[num2str(first),'-',num2str(last),'.tif'];
    sd=[save_dir,name];
    
    [ data ] = Tiff2d_2_data(files);
    ntheta=size(data,3);
    
else
    
    first = num2str(extract_number_from_string(files{1}));
    last = num2str(extract_number_from_string(files{2}));

    %if iscell(first) ~= 1,first={first};end
    %if iscell(last) ~= 1,last={last};end

    num_length=length(last);    %length required when building up the file name

    file0=strtrim(char(files{1}));                   %get the prefix to generate other file names
    prefix=file0(1:length(file0)-length(char(last))-length('.tif'));

    last=str2num(char(last));       %convert it to a number
    first=str2num(char(first));
    total_files=last-first+1;       %number of files that will be loaded

    %determine any leading zeros ie first=0000 last=0033, two leading zeros
    zero_pre=make_zeros(num_length-length(num2str(last)));
    prefix=[prefix,zero_pre];

    temp = imread(file0, 'tif');

    ntheta=total_files/accum;

    data=zeros([size(temp),ntheta]);    %create data array
    nx=size(data);
    temp=0;
    counter=0;

    name=[num2str(first),'-',num2str(last),'.tif'];

    try
        sd=[save_dir,name];
    catch
        sd=name;
    end

    for qq = 1:ntheta

        temp0=0;

        for pp = 1:accum

            number=num2str(first+counter);

            number=[make_zeros(length(num2str(last))-length(number)),number];

            fname=[prefix,number,'.tif'];
            disp([fname,' [',num2str(round(100.0*qq/(ntheta*accum))),'%]'])

            temp = imread(fname, 'tif');

            temp0=temp0+temp;
            temp=0;
            counter=counter+1;
        end

        data(:,:,qq)=temp0;

    end 

end %end of istruct
    
%do any orientation change if necessary
if isempty(orient) ~= 1
    
    switch orient
        case 'APS-pil'
            orient_m=[-1,1,0,0];
    end          
    data=change_orientation(data,orient_m);
end


%save to 3d tiff
for qq = 1:ntheta
    
    temp0=data(:,:,qq);
    
    if qq == 1,
        imwrite2tif(uint32(temp0),[],sd,'uint32','w');
    else
        imwrite2tif(uint32(temp0),[],sd,'uint32','a');
    end
    
end

disp('-----------------------------------')
disp(['File written to - ',name])
disp('-----------------------------------')

end

function zs_string=make_zeros(zs)

zs_string=''; 
    
if zs > 0
    for qq=1:zs
        
        zs_string=[zs_string,'0'];

    end

end
    
    
end

function data = change_orientation(data,orient_m)
    
    if numel(orient_m) == 0,orient_m=[0 0 0 0];end
    
    if abs(orient_m(1)) == 1,data=rot3d(data,orient_m(1));end

    if orient_m(2) == 1,data=flipdim(data,2);end
        
    if orient_m(3) == 1,data=flipdim(data,1);end
        
    if orient_m(4) == 1,data=flipdim(data,3);end
    
end

function [first last im_dir] = DLS_extract_fnums(fname,prefix)
%extracts .tif file numbers from .dat scan file from Diamond (I-07)

%thing to look fir to get relative dir
im_prefix='/';

try
    prefix;
catch
    prefix='p100kImage';         %set a default prefix
end

num_format='000000';              %used to determine length of f names      
file_format='.tif';
disp('______________________________________________________________')
disp('Extracting file numbers from .dat DLS scan file....')
disp(['Looking for files of the form - ',prefix,num_format,file_format])

fid=fopen(fname);           %open file
string=fscanf(fid,'%c');    %import as string
fclose(fid);            %close


%get locale of images
locs=strfind(string,prefix);

first=string(locs(1):locs(1)+numel(prefix)+numel(file_format)+numel(num_format));
last=string(locs(end):locs(end)+numel(prefix)+numel(file_format)+numel(num_format));
disp('')
disp(['[First,Last] -> [',strtrim(first),',',strtrim(last),']'])
disp('______________________________________________________________')
disp('')

%try
%get locale of diamond dir
try
    locs_im_st=strfind(string,im_prefix);
    im_dir=string(locs_im_st(end-1)+1:locs_im_st(end)-1);   
catch
    disp('Error determining image dir, using default....')
    im_dir=[];
end




end



