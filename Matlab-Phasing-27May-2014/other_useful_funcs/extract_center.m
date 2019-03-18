function new = extract_center(array,nx,ny,nz) 


dims=ndims(array);

if dims == 2, 
    new=zeros(nx,ny);
    nnx=size(array,2);
    nny=size(array,1);
    i=round(nny/2);
    j=round(nnx/2);
end

if dims == 3, 
    new=zeros(nx,ny,nz);
    nnx=size(array,2);
    nny=size(array,1);
    nnz=size(array,3);
    i=round(nny/2);
    j=round(nnx/2);
    k=round(nnz/2);
end

if dims == 2

    if mod(nx,2) == 1,new=array(i-(ny-1)/2:i+(ny-1)/2,j-(nx-1)/2:j+(nx-1)/2);end
    if mod(nx,2) == 0,new=array(i-(ny)/2:i+(ny)/2-1,j-(nx)/2:j+(nx)/2-1);end
end

if dims == 3

    if mod(nx,2) == 1,new=array(i-(ny-1)/2:i+(ny-1)/2,j-(nx-1)/2:j+(nx-1)/2,k-(nz-1)/2:k+(nz-1)/2);end
    if mod(nx,2) == 0,new=array(i-(ny)/2:i+(ny)/2-1,j-(nx)/2:j+(nx)/2-1,k-(nz)/2:k+(nz)/2-1);end
end


end

