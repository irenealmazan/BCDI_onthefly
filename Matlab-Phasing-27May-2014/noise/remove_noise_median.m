function [ data ] = remove_noise_median(data,fsize)
%jclark
%removes shot noise 
%does this by doing a median filter then uses all values >0
%from the fileterd version as a mask.  so isolated values are removed.
try
    fsize;
catch
    fsize=[2 2];
end

dims=ndims(data);

if dims == 3
    
    if numel(fsize) == 2
        for qq=1:size(data,3)      
            mask=medfilt2(data(:,:,qq),fsize);
            mask(mask > 0)=1;
            data(:,:,qq)=data(:,:,qq).*mask;
        end
    end
    if numel(fsize) == 3
        mask=medfilt3(data,fsize);
        mask(mask > 0)=1;
        data=data.*mask;
    end
end
if dims == 2
   
   mask=medfilt2(data,fsize);
   data(mask == 0)=0; 
   
end



end

