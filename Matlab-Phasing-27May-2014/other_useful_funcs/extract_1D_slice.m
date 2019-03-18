function [ slice ] = extract_1D_slice(array,orient,n1,n2 )
%extracts a slice from a 3d volume. orient is the
%slice to be extracted, i.e 'xy','xz' et etc
%nc is the slice number.
%returns 2d slice.  normal image x-y is defined as x->number of columns and
%y number of rows.  return slice is this convention where columns is given
%by first letter and rows second


switch lower(orient)
    
    case 'x'
        
        slice=squeeze(array(n1,:,n2));
        slice=slice';
    case 'y'
        
        slice=squeeze(array(:,n1,n2));
        
    case 'z' 
        
        slice=(squeeze(array(n2,n1,:)));
        
end
        
end

