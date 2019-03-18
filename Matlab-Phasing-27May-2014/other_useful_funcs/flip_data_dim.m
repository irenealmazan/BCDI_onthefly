function [ data ] = flip_data_dim(data,flips)
%jclark

if flips(1) == 1,data=flipdim(data,1);end
if flips(2) == 1,data=flipdim(data,2);end
if flips(3) == 1,data=flipdim(data,3);end

end

