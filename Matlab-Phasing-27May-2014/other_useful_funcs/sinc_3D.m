function fxy = sinc_3D(nx,ny,nz,fwhmx,fwhmy,fwhmz)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians


[x ,y,z]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2,-(nz-1)/2:(nz-1)/2);

x=x/max(x(:));
y=y/max(y(:));
z=z/max(z(:));

x=x*.6047/fwhmx;
y=y*.6047/fwhmy;
z=z*.6047/fwhmz;

gx=sinc(x);
gy=sinc(y);
gz=sinc(z);

fxy=gx.*gy.*gz;

fxy=fxy/sum(fxy(:));

end

