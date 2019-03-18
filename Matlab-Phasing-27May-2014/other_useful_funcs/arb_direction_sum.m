function [fx fy rx] = arb_direction_sum(data,theta)
%sum's along direction theta and its orthogonal direction
use_imrotate=1;

[sy sx]=size(data);

%create the orthogonal directions
[x y]=meshgrid(1:sx,1:sy);

%create the rotated directions
u=x*cos(deg2rad(theta))-sin(deg2rad(theta))*y;
v=x*sin(deg2rad(theta))+cos(deg2rad(theta))*y;

fx=0;
fy=0;
ry=zeros(size(data));
rx=zeros(size(data));
if use_imrotate ~= 1
    for qq=1:sy

       xc=round(u(qq,:));
       yc=round(v(qq,:));

       for tt=1:numel(xc)
          try
              yline(tt)=data(yc(tt),xc(tt));
          catch
              yline(tt)=0;
          end
       end

       rx(qq,:)=yline;
       fx=fx+yline;

    end

    for qq=1:sx

       xc=round(u(:,qq));
       yc=round(v(:,qq));

       for tt=1:numel(yc)
          try
              xline(tt)=data(yc(tt),xc(tt));
          catch
              xline(tt)=0;
          end
       end

       ry(:,qq)=xline;
       fy=fy+xline;


    end
else
    rx=imrotate(data,theta,'nearest','crop');
    fx=sum(rx,1);
    fy=sum(rx,2);
    
end

end

