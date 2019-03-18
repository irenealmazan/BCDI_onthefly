function [ new_array ] = rot3darb( array,k,mth )
%jclark
%rotates a 3d array around the x-y axis
%slow but i couldn't be bothered making it any faster
%k is integer amounts of rotation counter-clockwise

if exist('k') == 0, k = 0;end

nn=size(array);

if exist('mth') == 0,mth='nearest';end %'nearest'

if numel(nn) == 3
    
    new_array=zeros([nn(2),nn(1),nn(3)]);

    for ww=1:nn(3), new_array(:,:,ww)=imrotate(squeeze(array(:,:,ww)),k,mth,'crop');end

end

if numel(nn) == 2
    
    new_array=imrotate(array,k,mth,'crop');

end




end

