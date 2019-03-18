function [ xyz] = center_of_mass_v2( array )
%jclark

array=double(array);

tot=sum(sum(sum(array)));

sz=size(array);


nx=sz(2);
ny=sz(1);
if ndims(array) == 3,nz=sz(3);end

xyz=0;

if ndims(array) == 3
    [x , y,z]=meshgrid( 1:nx,1:ny,1:nz);
end
if ndims(array) == 2
    [x , y]=meshgrid( 1:nx,1:ny);
end

xyz(1)=sum(sum(sum(array.*x)))/tot;
xyz(2)=sum(sum(sum(array.*y)))/tot;

if ndims(array) ==3, xyz(3) = sum(sum(sum(array.*z)))/tot;

end


