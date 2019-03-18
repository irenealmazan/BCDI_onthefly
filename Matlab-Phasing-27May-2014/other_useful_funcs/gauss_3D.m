function fxy = gauss_3D(nx,ny,nz,sigx,sigy,sigz,alpha)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians
try
    alpha;
catch
    alpha=1;
end

fxy=zeros(nx,ny,nz);


[x , y,z]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2,-(nz-1)/2:(nz-1)/2);

% if angle ~= 0
%     theta=angle/180e0*pi;
%     xd=x*cos(theta)-y*sin(theta);
%     yd=x*sin(theta)+y*cos(theta);
%     x=xd;
%     y=yd;
%     xd=0;
%     yd=0;
% end

if mod(nx,2) == 1
    if sigx == 0, gx=(x == 0);else gx=exp(-alpha*0.5.*x.^2./sigx^2);end
else gx=exp(-alpha*0.5.*x.^2./sigx^2);end

if mod(ny,2) == 1
    if sigy == 0, gy=(y == 0);else gy=exp(-alpha*0.5.*y.^2./sigy^2);end
else gy=exp(-alpha*0.5.*y.^2./sigy^2);end

if mod(nz,2) == 1
    if sigz == 0, gz=(z == 0);else gz=exp(-alpha*0.5.*z.^2./sigz^2);end
else gz=exp(-alpha*0.5.*z.^2./sigz^2);end


fxy=gx.*gy.*gz;

fxy=fxy/sum(((fxy(:))));

end


