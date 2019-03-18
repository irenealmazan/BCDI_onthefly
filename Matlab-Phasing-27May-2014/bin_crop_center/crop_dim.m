function [ cropped ] = crop_dim(array,new,dim )
%jclark

sz=size(array);

%switch dim
if dim == 1
    cropped=array(:,new(1):new(2),:);end
if dim == 2 
    cropped=array(new(1):new(2),:,:);end
if dim == 3
    cropped=array(:,:,new(1):new(2));end
%end

end

