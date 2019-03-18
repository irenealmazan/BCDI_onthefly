function [ slice ] = extract_3D_slice(array,orient,nc )
%extracts a slice from a 3d volume. orient is the
%slice to be extracted, i.e 'xy','xz' et etc
%nc is the slice number.
%returns 2d slice.  normal image x-y is defined as x->number of columns and
%y number of rows.  return slice is this convention where columns is given
%by first letter and rows second

    
switch lower(orient)
    
    case 'xy'
        try
            nc;
        catch
            nc=ceil(size(array,3)/2);
        end
        slice=squeeze(array(:,:,nc));
    
    case 'zx'
        try
            nc;
        catch
            nc=ceil(size(array,1)/2);
        end
        slice=squeeze(array(nc,:,:));
        
    case 'xz' 
        try
            nc;
        catch
            nc=ceil(size(array,1)/2);
        end
        slice=rot90(squeeze(array(nc,:,:)));
        
    case 'yx'
        try
            nc;
        catch
            nc=ceil(size(array,3)/2);
        end
        slice=rot90(squeeze(array(:,:,nc)));

    case 'zy'
        try
            nc;
        catch
            nc=ceil(size(array,2)/2);
        end
        slice=squeeze(array(:,nc,:));
        
    case 'yz'
        try
            nc;
        catch
            nc=ceil(size(array,2)/2);
        end
        slice=rot90(squeeze(array(:,nc,:)));
        
end
        
end

