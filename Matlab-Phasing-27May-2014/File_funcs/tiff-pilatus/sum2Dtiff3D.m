function sum2Dtiff3D( files,accum,save_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%this takes a series of 2d tiff images and outputs to a .mat file
%specify the first and last files 
%use accum to specify the accumulations
%to save as a .tif stack use Tiff2D_2_Tiff3D

idx = regexp(char(files{1}),'\d+');
nums = regexp(char(files{1}),'\d+','match');   %get the numbers from the file name

first=char(nums(end));        %get first 

idx = regexp(char(files{2}),'\d+');
nums = regexp(char(files{2}),'\d+','match');   %get the numbers from the file name

last=char(nums(end));

file0=char(files{1});                   %get the prefix to generate other file names
prefix=file0(1:length(file0)-length(char(nums(end)))-length('.tif'));

last=str2num(char(last));       %convert it to a number
first=str2num(char(first));
total_files=last-first+1;       %number of files that will be loaded


temp = imread(file0, 'tif');

ntheta=total_files/accum;

data=zeros([size(temp),ntheta]);    %create data array
nx=size(data);
temp=0;
counter=0;

for qq = 1:ntheta

    temp0=0;

    for pp = 1:accum

        fname=[prefix,num2str(first+counter),'.tif'];
        disp(fname(end-10:end))

        temp = imread(fname, 'tif');

        temp0=temp0+temp;
        temp=0;
        counter=counter+1;
    end

    data(:,:,qq)=temp0;

end

name=[num2str(first),'-',num2str(last),'.mat'];

try
    save([save_dir,name],'data');
catch
    save(name,'data');
end

end

