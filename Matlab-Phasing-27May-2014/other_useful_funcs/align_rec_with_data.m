function [pn] = align_rec_with_data(pn,data)
%jclark
%will align a reconstruction with a data set
%correct pixel number etc.

if ischar(pn) == 1
   
    pn=load_rec_from_dir([pn],1);
    
end

pn=zero_pad_ver3(pn,size(data,2),size(data,1),size(data,3));  %equal pixel numbers/binning

F1=ifftshift(fftn(fftshift(pn)));  %estimate data

F1=F1*sqrt(sum(abs(data(:).^2))/sum(abs(F1(:)).^2));  %equal power

disp('Aligning with data....')

[h k l]=register_3d_reconstruction(abs(data).^2,abs(F1).^2); %align with the data

pn=ifftshift(ifftn(fftshift(sub_pixel_shift(F1,h,k,l))));


end

