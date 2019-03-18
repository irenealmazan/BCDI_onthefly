function [ output_args ] = mat2tif(data,save_name)
%jclark

disp('Saving as .tif file....')
    
for qq=1:size(data,3)
    temp0=data(:,:,qq);
    if qq == 1,
        imwrite2tif(uint32(temp0),[],save_name,'uint32','w');
    else
        imwrite2tif(uint32(temp0),[],save_name,'uint32','a');
    end
    clear temp0
end


end

