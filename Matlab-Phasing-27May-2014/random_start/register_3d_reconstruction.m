function [h k l] = register_3d_reconstruction(a,b)
%jclark
%returns the hkl required to shift b to a
%e.g a=circshift(b,[h,k,l])

hk=dftregistration(fft2(squeeze(sum(a,3))),fft2(squeeze(sum(b,3))),100);
hl=dftregistration(fft2(squeeze(sum(a,2))),fft2(squeeze(sum(b,2))),100);
kl=dftregistration(fft2(squeeze(sum(a,1))),fft2(squeeze(sum(b,1))),100);

if ndims(a) > 2
    h=sum([hk(3),hl(3)])*0.5;
    k=sum([hk(4),kl(3)])*0.5;
    l=sum([hl(4),kl(4)])*0.5;
else
    h=hk(3);
    k=hk(4);
    l=0;
end
    
    
end

