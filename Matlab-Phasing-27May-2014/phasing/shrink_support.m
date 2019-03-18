function [shrink fill] = shrink_support(support,shrink )
%jclark


com_supp=center_of_mass(support);   %get center of mass to align other shrunk copies



shrinkage=shrink;     %get the shrink increments
scale=shrinkage;

nn=size(support);

nx=nn(2);
ny=nn(1);
if numel(nn) == 3,nz=nn(3);end

fill=support;


if numel(nn) == 2
    shrunk=zero_pad_ver2( imresize(support,scale,'nearest'),nx,ny);
end
if numel(nn) == 3
    shrunk=zero_pad_ver2( imresize(support,scale,'nearest'),nx,ny,nz);
end

com_shrunk=center_of_mass(round(abs(shrunk)/max(abs(shrunk(:)))));
delta_com=round(com_supp-com_shrunk);

if numel(delta_com) ==3,shrunk=circshift(shrunk,[delta_com(2),delta_com(1),delta_com(3)]);end
if numel(delta_com) ==2,shrunk=circshift(shrunk,[delta_com(2),delta_com(1)]);end

fill=abs(shrunk);


shrink=round(fill/max(fill(:)));


end

