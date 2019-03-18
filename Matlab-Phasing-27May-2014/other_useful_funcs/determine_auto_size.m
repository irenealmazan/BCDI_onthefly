function [x y z] = determine_auto_size(auto)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dims=ndims(auto);

sz=size(auto);

xc=round(sz(2)/2);
yc=round(sz(1)/2);

fact=.2;

if dims == 3
    
    %%
    zc=round(sz(3)/2);

    xz=extract_3D_slice(abs(auto),'xz');
    
    [fx fz]=gradient(abs(xz));
   
    xx=abs(fx(zc,:));
    zz=abs(fz(:,xc));

    indx=find(xx > fact*max(xx(:)));
    indz=find(zz > fact*max(zz(:)));

    xl1=indx(end)-indx(1);
    zl1=indz(end)-indz(1);
    indx=[];
    indz=[];
    fx=[];
    fz=[];
    
    xx=[];
    zz=[];
    %%   
    zy=extract_3D_slice(abs(auto),'zy');
    [fz fy]=gradient(abs(zy));
    zz=abs(fz(zc,:));
    yy=abs(fy(:,yc));

    indz=find(zz > fact*max(zz(:)));
    indy=find(yy > fact*max(yy(:)));

    zl2=indz(end)-indz(1);
    yl1=indy(end)-indy(1);
    indz=[];
    indy=[];
    fz=[];
    fy=[];
    
    yy=[];
    zz=[];
end
    %%   
    if dims == 3,xy=extract_3D_slice(abs(auto),'xy');else xy=abs(auto);end
    [fx fy]=gradient(abs(xy));
    xx=abs(fx(xc,:));
    yy=abs(fy(:,yc));

    indx=find(xx > fact*max(xx(:)));
    indy=find(yy > fact*max(yy(:)));

    xl2=indx(end)-indx(1);
    yl2=indy(end)-indy(1);
    indx=[];
    indy=[];
    fx=[];
    fy=[];
    xy=[];
    yy=[];
    xx=[];
    zy=[];
    xz=[];
% else
%     
%    xy=auto;
%    [fx fy]=gradient(abs(xy));
%    
%     xx=fx(yc,:);
%     yy=fy(:,xc);
% 
%     indx=find(xx > fact*max(xx(:)));
%     indy=find(yy > fact*max(yy(:)));
% 
%     xl=indx(end)-indx(1);
%     yl=indy(end)-indy(1);
   




end

