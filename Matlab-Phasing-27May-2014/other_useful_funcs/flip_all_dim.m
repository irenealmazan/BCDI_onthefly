function [ array ] = flip_all_dim(array)
%jclark

array=flipdim(array,1);

array=flipdim(array,2);

if ndims(array) == 3, array=flipdim(array,3);end



end

