function [shifted]=sub_pixel_shift(array,row_shift,col_shift,z_shift)
%jclark adapted from M.G-S

if ndims(array) == 3
    buf2ft=fftn(array);
    [nr,nc,nz]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
    [Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
    %Greg = Greg*exp(i*diffphase);
    shifted=ifftn(Greg);
else
    buf2ft=fftn(array);
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(1i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    %Greg = Greg*exp(i*diffphase);
    shifted=ifftn(Greg);  
end


end
