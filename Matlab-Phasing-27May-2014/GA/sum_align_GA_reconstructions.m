function  [summed sigma pn_g]=sum_align_GA_reconstructions(pn_g,ind)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

try
    ind;
catch
    ind=1;
end

pn_g=align_iterates(pn_g,ind,0);

sz=size(pn_g);

%summed=mean(pn_g,ndims(pn_g));
%sigma=std(pn_g,1,ndims(pn_g));
summed=mean(real(pn_g),ndims(pn_g))+i*mean(imag(pn_g),ndims(pn_g));
sigma=std(real(pn_g),1,ndims(pn_g))+i*std(imag(pn_g),1,ndims(pn_g));


    

end

