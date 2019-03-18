%% this is just to create a 2d gauss for smoothing the square sample
function fxy = super_gauss_2D(nx,ny,sigx,sigy,pow )
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians
try
    pow;
catch
    pow=1;
end


if isempty(sigy), sigy = sigx;end

fxy=zeros(nx,ny);


[x , y]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2);

if sigx ~= sigy
    if mod(nx,2) == 1
        if sigx == 0, gx=(x == 0);else gx=exp(-0.5.*(x.^2./sigx^2).^pow);end
    else gx=exp(-0.5.*(x.^2./sigx^2).^pow);end

    if mod(ny,2) == 1
        if sigy == 0, gy=(y == 0);else gy=exp(-0.5.*(y.^2./sigy^2).^pow);end
    else gy=exp(-0.5.*(y.^2./sigy^2).^pow);end

    fxy=gx.*gy;
else
    
    rr=x.^2+y.^2;
    
    fxy=exp(-0.5.*(rr/sigx^2).^pow);
    
end


fxy=fxy/sum(sum(fxy));

end
