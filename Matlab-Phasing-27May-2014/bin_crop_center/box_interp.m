function [ barray ] = box_interp( array,bx,by,med )
%jclark

try
    med;
catch
    med=0;
end

sz=size(array);

x0=sz(2);
y0=sz(1);


while mod(x0,bx) ~=0,x0=x0+1;end
while mod(y0,by) ~=0,y0=y0+1;end


array=padarray(array,[(y0-sz(1)),(x0-sz(2)),0],0,'pre');

newx=x0/bx;
newy=y0/by;

barray=zeros([by,bx]);

x00=[1:bx:x0,x0+1];
y00=[1:by:y0,y0+1];

if bx == 1 && by == 1
    barray=array;
    %disp(['No binning since bx = ',num2str(bx),' & by = ',num2str(by)])
else
    for xx = 1:newx
        for yy = 1:newy

            if med == 1
                temp=(( array(y00(yy):y00(yy+1)-1, x00(xx):x00(xx+1)-1  )));
                barray(yy,xx)=median(temp(:));
            else
                barray(yy,xx)=sum(sum( array(y00(yy):y00(yy+1)-1, x00(xx):x00(xx+1)-1  )));
            end
        end
    end
end

end

