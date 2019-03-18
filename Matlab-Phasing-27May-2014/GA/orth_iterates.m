function [ pn_g ] = orth_iterates(pn_g)
%jclark

tot=sum(abs(pn_g(:)));

pn_g=align_iterates(pn_g);


north_ill=reshape(pn_g,[prod(size(pn_g(:,:,1))),size(pn_g,3)]);

north_ill=orth(north_ill);

pn_g=reshape(north_ill,[size(pn_g)]);

pn_g=pn_g/sum(abs(pn_g(:)))*tot;

end

