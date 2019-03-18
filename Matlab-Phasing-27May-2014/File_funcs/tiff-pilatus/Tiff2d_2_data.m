function [ data ] = Tiff2d_2_data(files)
%jclark
%loads 2d tiff files into 3d matlab array
%will load all files that are passed through
%must use rdir to find files recursively then pass through this result

nfiles=max(size(files));  %get number of files

for qq=1:nfiles
    
   disp( ['[',num2str(round(100.0*qq/(nfiles))),'%]'])
   
   temp=imread(char(files(qq).name),'tif'); 
    
   data(:,:,qq)=temp;
   
   temp=[];
    
end



end

