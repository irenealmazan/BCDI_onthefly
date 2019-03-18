function [ vol ] = make_slice_3D(slice,vol,orient,nc )
%BE CAREFUL IF THE SLICES ARE NOT SYMMETRIC
%nc is where it should be placed
%inserts the slice into a 3d volume
%orient is used to provide the orientation of the slice
%vol is the volume to inert into


switch lower(orient)
    
    case 'xy'
        vol(:,:,nc)=slice;
        %slice=squeeze(array(:,:,nc));
    
    case 'zx'
        vol(nc,:,:)=slice;
        %slice=squeeze(array(nc,:,:));
        
    case 'xz' 
        vol(nc,:,:)=rot90(slice,-1);
        %slice=rot90(squeeze(array(nc,:,:)));
        
    case 'yx'
        vol(:,:,nc)=rot90(slice,-1);
        %slice=rot90(squeeze(array(:,:,nc)));

    case 'zy'
        vol(:,nc,:)=slice;
        %slice=squeeze(array(:,nc,:));
        
    case 'yz'
        vol(:,nc,:)=rot90(slice,-1);
        %slice=rot90(squeeze(array(:,nc,:)));
        
end
        
end