function [ data ] = loadTif3D(file)
%jclark
%loads a multipage tiff

temp = imfinfo(file, 'tif');
            
num_images=numel(temp);
for qq = 1:num_images

    ff0 = imread(file,qq);

    data(:,:,qq) = ff0;

end


end
