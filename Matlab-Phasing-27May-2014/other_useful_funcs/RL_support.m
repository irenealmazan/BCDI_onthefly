function [ rk] = RL_support( Im,Ik,psf)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Imk=ifftshift(ifftn(fftshift(abs(fftshift(fftn(ifftshift(Ik)))).^2)));%convnfft(Ik,psf,'same');

dk=zeros(size(psf));

dk(psf ~=0)=Im(psf ~=0)./Imk(psf ~=0);

H0=sum(abs(Ik(:)));

rk=convnfft(conj_reflect(Ik),dk,'same')/H0;

end

