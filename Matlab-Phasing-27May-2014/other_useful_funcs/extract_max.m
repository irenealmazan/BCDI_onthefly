function [ new ] = extract_max( array,nx,ny,nz )
%cuts out an array of size nxmny centered arond the max value
%   Detailed explanation goes here
if ndims(array) == 2
    thing=sum(size(array)-[ny,nx]);
end

if ndims(array) == 3
    thing=sum(size(array)-[ny,nx,nz]);
end


if thing ~= 0

    array=center_array(array); %jnc feb2011

    dims=ndims(array);

    if dims == 2, 
        new=zeros(nx,ny);
        [i,j]=ind2sub(size(array),find( abs(array) == max(max(max(abs(array))))));
    end

    if dims == 3, 
        new=zeros(nx,ny,nz);
        [i,j,k]=ind2sub(size(array),find( abs(array) == max(max(max(abs(array))))));
    end

    xyz=center_of_mass(abs(array));     %calc com 
    sz=size(array);

    if dims ==3
        %while sum(abs( (xyz+0.5*([sz(2),sz(1),sz(3)]))-[j,i,k])) >= nx/2, 
         %   disp('WARNING - com and max val disagree....')
          %  disp('comparing next highest value...')
           % disp('recalculating....')
           % array(i,j,k)=0;
           % [i,j,k]=ind2sub(size(array),find( abs(array) == max(max(max(abs(array))))));
           % xyz=center_of_mass(abs(array));     %calc com 
        %end
    end

    if dims == 2

        if mod(nx,2) == 1,new=array(i-(ny-1)/2:i+(ny-1)/2,j-(nx-1)/2:j+(nx-1)/2);end
        if mod(nx,2) == 0,new=array(i-(ny)/2:i+(ny)/2-1,j-(nx)/2:j+(nx)/2-1);end
    end

    if dims == 3

        if mod(nx,2) == 1,new=array(i-(ny-1)/2:i+(ny-1)/2,j-(nx-1)/2:j+(nx-1)/2,k-(nz-1)/2:k+(nz-1)/2);end
        if mod(nx,2) == 0,new=array(i-(ny)/2:i+(ny)/2-1,j-(nx)/2:j+(nx)/2-1,k-(nz)/2:k+(nz)/2-1);end
    end

else
   new=array; 
end

end

