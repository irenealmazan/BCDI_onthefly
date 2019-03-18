function edf2D_2_Tiff3D( files,accum,save_dir,orient)
%takes a series of 2D tiff's and saves as a 3d tiff.
%use sum2Dtiff3D.m to save as .mat file.
%specify the first and last files eg. files={first_file,last_file}
%use accum to specify the accumulations
%orient is a 4 vector orientation thing

if exist('orient') ~= 1,orient=[1 0 1 0];end

ftype='.edf';                   %specify the suffix
%%
%get file numbers, assumes that it is contiguous
idx = regexp(char(files{1}),'\d+');
nums = regexp(char(files{1}),'\d+','match');   %get the numbers from the file name

first=strtrim(char(nums(end)));        %get first 

%get the last number
idx = regexp(char(files{2}),'\d+');
nums = regexp(char(files{2}),'\d+','match');   %get the numbers from the file name

last=strtrim(char(nums(end)));

file0=char(files{1});                   %get the prefix to generate other file names
prefix=file0(1:length(file0)-length(char(nums(end)))-length(ftype));

last=str2num(char(last));       %convert it to a number
first=str2num(char(first));
total_files=last-first+1;       %number of files that will be loaded
%%

zero_pre='';
nlast=numel(num2str(last));

while nlast ~= numel(char(nums(end)))
    zero_pre=[zero_pre,'0'];
    nlast=nlast+1;
end

prefix=[prefix,zero_pre];

%temp = imread(file0, 'tif');
[header temp]=pmedf_read(file0);

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

if isdir(save_dir) ~= 1,mkdir(save_dir);end

for qq = 1:ntheta

    temp0=0;

    for pp = 1:accum

        number=num2str(first+counter);

        number=[make_zeros(length(num2str(last))-length(number)),number];

        
        %;fname=[prefix,num2str(first+counter),ftype];
        fname=[prefix,number,ftype];
        
        disp([fname,' [',num2str(round(100.0*qq/(ntheta*accum))),'%]'])

        %temp = imread(fname, 'tif');
        [header1 temp]=pmedf_read(fname);
        
        temp0=temp0+temp;
        temp=0;
        counter=counter+1;
        clear header1
    end

    data(:,:,qq)=temp0;
    
    
%     if qq == 1,
%         imwrite2tif(uint32(temp0),[],sd,'uint32','w');
%     else
%         imwrite2tif(uint32(temp0),[],sd,'uint32','a');
%     end
end 

%do any orientation change if necessary
if isempty(orient) ~= 1
            
    data=change_orientation(data,orient);
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

