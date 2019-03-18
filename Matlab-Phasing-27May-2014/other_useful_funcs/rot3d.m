function [ new_array ] = rot3d( array,k )
%jclark
%rotates a 3d array around the x-y axis
%slow but i couldn't be bothered making it any faster
%k is integer amounts of rotation counter-clockwise

if exist('k') == 0, k = 1;end

nn=size(array);


if numel(nn) == 3
    
    new_array=zeros([nn(2),nn(1),nn(3)]);

    for ww=1:nn(3), new_array(:,:,ww)=rot90(squeeze(array(:,:,ww)),k);end

end

if numel(nn) == 2
    
    new_array=rot90(array,k);

end




end

