function [ pn0] = remove_ramp_pn(pn)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

data=(fftshift(fftn(ifftshift(pn))));

[hkl]=cent_arr_auto(abs(data).^2,1);

h=hkl(1);
k=hkl(2);
l=hkl(3);

pn0=fftshift(ifftn(fftshift(circshift(data,[-h,-k,-l]))));


% data=(fftshift(fftn(ifftshift(pn0))));
% 
% [hkl]=cent_arr_auto(abs(data).^2,2);
% 
% h=hkl(1);
% k=hkl(2);
% l=hkl(3);
% 
% data=sub_pixel_shift(data,-h,-k,-l);
% 
% pn0=fftshift(ifftn(fftshift((data))));




end

