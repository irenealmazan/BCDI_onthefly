function [ fill ] = fill_support(support,steps )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

try
    steps;
catch
    steps=9;
end

com_supp=center_of_mass(support);   %get center of mass to align other shrunk copies

min_shrink=0.9;                 %set the min shrink, prevents the support 
                                %from growing
max_shrink=0.1;                 %set the most it will shrink by.  as a fraction
                                
shrink_inc=(min_shrink-max_shrink)/steps;

shrinkage=max_shrink:shrink_inc:min_shrink;     %get the shrink increments

nn=size(support);

nx=nn(2);
ny=nn(1);
if numel(nn) == 3,nz=nn(3);end

fill=support;

for qq=1:steps
    
    scale=shrinkage(qq);
    
    if numel(nn) == 2
        shrunk=zero_pad_ver2( imresize(support,scale),nx,ny);
    end
    if numel(nn) == 3
        shrunk=zero_pad_ver2( imresize(support,scale),nx,ny,nz);
    end
    
    com_shrunk=center_of_mass(ceil(abs(shrunk)/max(abs(shrunk(:)))));
    delta_com=round(com_supp-com_shrunk);
    
    if numel(delta_com) ==3,shrunk=circshift(shrunk,[delta_com(2),delta_com(1),delta_com(3)]);end
    if numel(delta_com) ==2,shrunk=circshift(shrunk,[delta_com(2),delta_com(1)]);end
    
    fill=fill+abs(shrunk);
    if qq == 1,disp('Filling support....');end
end

fill=ceil(fill/max(fill(:)));


end

