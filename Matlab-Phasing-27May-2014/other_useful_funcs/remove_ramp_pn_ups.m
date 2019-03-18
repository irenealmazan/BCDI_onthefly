function [ pn0] = remove_ramp_pn_ups(pn,upsample)
%jclark
if exist('upsample') ~= 1,upsample=5;end

if ndims(pn) == 3
    new_s=upsample*[size(pn,2),size(pn,1),size(pn,3)];
     pn=zero_pad_ver_Irene(new_s(1),new_s(2),new_s(3),pn);
    %pn=zero_pad_ver3(pn,new_s(1),new_s(2),new_s(3));
else
    new_s=upsample*[size(pn,2),size(pn,1)];
    pn=zero_pad_ver3(pn,new_s(1),new_s(2));
end



data=(fftshift(fftn(ifftshift(pn))));

%%
xyz=center_of_mass(abs(data).^4);
xyz=xyz-.5;
%%
% [ array yxz] = center_array(abs(data).^2);
% yxz=-1*(yxz+1);
% xyz=yxz;
% xyz(1)=yxz(2);
% xyz(2)=yxz(1);
%%

h=xyz(2);
k=xyz(1);

if ndims(pn) == 3
    l=xyz(3);
    pn0=fftshift(ifftn(fftshift(sub_pixel_shift(data,-h,-k,-l))));
else
    pn0=fftshift(ifftn(fftshift(sub_pixel_shift(data,-h,-k))));    
end

if ndims(pn) == 3
    pn0=zero_pad_ver3(pn0,new_s(1)/upsample,new_s(2)/upsample,new_s(3)/upsample);
else
    pn0=zero_pad_ver3(pn0,new_s(1)/upsample,new_s(2)/upsample);
end

end