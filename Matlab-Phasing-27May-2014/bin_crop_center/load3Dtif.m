function [ data] = load3Dtif( fname )
%loads a 3d tiff into array 

temp = imfinfo(fname, 'tif');
            
num_images=numel(temp);
for qq = 1:num_images

    ff0 = imread(fname,qq);

    if qq == 1, data=zeros([size(ff0),num_images]);end

    data(:,:,qq) = ff0;

end



end

