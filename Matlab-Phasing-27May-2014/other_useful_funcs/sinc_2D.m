function fxy = sinc_2D(nx,ny,fwhmx,fwhmy)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians


[x , y]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2);

x=x/max(x(:));
y=y/max(y(:));

x=x*.6047/fwhmx;
y=y*.6047/fwhmy;

gx=sinc(x);
gy=sinc(y);

fxy=gx.*gy;

fxy=fxy/sum(sum(fxy));

end

