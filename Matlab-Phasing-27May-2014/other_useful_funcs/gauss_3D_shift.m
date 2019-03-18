function fxy = gauss_3D_shift(nx,ny,nz,sigx,sigy,sigz,mux,muy,muz)
%jclark

fxy=zeros([ny,nx,nz]);

[x , y,z]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2,-(nz-1)/2:(nz-1)/2);

x=x-mux;
y=y-muy;
z=z-muz;

fxy=exp(-0.5*(x/sigx).^2-0.5*(y/sigy).^2-0.5*(z/sigz).^2);

fxy=fxy/sum(((fxy(:))));

end



