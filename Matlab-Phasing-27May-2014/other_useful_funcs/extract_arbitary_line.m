function [ fxy ] = extract_arbitary_line(data,x1,y1,x2,y2)
%extract a line from a 2d image with starting points x1,y1 and end points
%x2,y2

if numel(x1) == 4
   y1=x1(2);
   x2=x1(3);
   y2=x1(4);
   x1=x1(1);
end


nx=x2-x1;
ny=y2-y1;

nn=max([abs(nx),abs(ny)]);

x=(1:nn)/nn;
x=x*nx+x1;

y=(1:nn)/nn;
y=y*ny+y1;

x=round(x);
y=round(y);

for qq=1:nn
    fxy(qq)=data(y(qq),x(qq));
end

end

