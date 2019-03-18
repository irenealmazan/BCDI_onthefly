function [lattice] = create_nd_lattice(ngpr,data)
%jclark
%creates a 2d or 3d latttice

if ndims(data) == 3
    [xx yy zz]=meshgrid(1:size(data,2),1:size(data,1),1:size(data,3)); %make attice for gaussians    
    
    ngprz=round(size(data,3)/ngpr);
    zz=mod(zz,ngprz);
    indz=(zz == round(ngprz/2));
end

if ndims(data) == 2 
    [xx yy]=meshgrid(1:size(data,2),1:size(data,1)); %make attice for gaussians 
    indz=1;
    
end

ngprx=round(size(data,2)/ngpr);
ngpry=round(size(data,1)/ngpr);


xx=mod(xx,ngprx);yy=mod(yy,ngpry);
indx=(xx == round(ngprx/2));indy=(yy == round(ngpry/2));

lattice=indx.*indy.*indz;

end

