function [ sym_kern ] = symmetrize_kernal( kernal )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

tot=sum(abs(kernal(:)));

%sym_kern=0.5*(kernal+conj_reflect(kernal));

sym_kern=kernal+flipdim(flipdim(flipdim(kernal,1),2),3);

sym_kern=sym_kern/sum(sym_kern(:))*tot;

end

