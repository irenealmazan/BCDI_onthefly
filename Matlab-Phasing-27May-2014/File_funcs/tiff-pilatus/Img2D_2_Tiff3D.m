function [ data ] = Img2D_2_Tiff3D( files,accum,save_dir)
%function Img2D_2_Tiff3D( files,accum,save_dir) old usage jnc july2012
%takes a series of 2D tiff's and saves as a 3d tiff.
%use sum2Dtiff3D.m to save as .mat file.
%specify the first and last files eg. files={first_file,last_file}
%use accum to specify the accumulations


first = num2str(extract_number_from_string(files{1}));
last = num2str(extract_number_from_string(files{2}));

%if iscell(first) ~= 1,first={first};end
%if iscell(last) ~= 1,last={last};end

num_length=length(last);    %length required when building up the file name

file0=char(files{1});                   %get the prefix to generate other file names
prefix=file0(1:length(file0)-length(char(last))-length('.tif'));

last=str2num(char(last));       %convert it to a number
first=str2num(char(first));
total_files=last-first+1;       %number of files that will be loaded

%determine any leading zeros ie first=0000 last=0033, two leading zeros
zero_pre=make_zeros(num_length-length(num2str(last)));
prefix=[prefix,zero_pre];

temp = imageread(file0, 'img',[487,195],32);

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
        
        fname=[prefix,number,'.img'];
        disp([fname,' [',num2str(round(100.0*qq/(ntheta*accum))),'%]'])

        temp = imageread(fname, 'img',[487,195],32);

        temp0=temp0+temp;
        temp=0;
        counter=counter+1;
    end

    data(:,:,qq)=temp0;
    
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


