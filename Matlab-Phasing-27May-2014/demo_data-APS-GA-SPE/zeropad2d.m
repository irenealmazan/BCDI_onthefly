function out = zeropad2d(in, outcol, outrow, coloffset, rowoffset)

%SH 8-3-09
%
%this function pads a two-d array with zeros, while the original array
%remains centered.  this is not how fft2 handles zero padding, which works
%by simply adding zeros in +x and +y only.

if nargin<4
    coloffset=0;
    rowoffset=0;
end

if (size(in, 2) > outcol)
    flag=1;
    while flag
        in(:,1)=[];
        if size(in,2) <= outcol break; end
        in(:,end) = [];
        if size(in,2) <= outcol break; end
%         if size(in,2) <= outcol
%             flag=0; 
%             display('shrink number of columns');
%         end
    end
end

if (size(in, 1) > outrow)
    flag=1;
    while flag
        in(1,:)=[];
        if size(in,1) <=outrow break; end
        in(end,:) = [];
        if size(in,1) <=outrow break; end
%         if size(in,1) <= outrow
%             flag=0; 
%             display('shrink number of rows');
%         end
    end
end

[inrow incol] = size(in);

if mod(inrow,2)~=0   
    inrow=inrow+1;  
    in(inrow,:)=zeros;
end

if mod(incol,2) ~=0 
    incol=incol+1; 
    in(:,incol)=zeros;
end

out(outrow, outcol)=0;

lowcol = outcol/2 - incol/2 +1;
highcol = outcol/2 + incol/2;

lowrow = outrow/2 - inrow/2+1;
highrow = outrow/2 + inrow/2;

out([lowrow:highrow] + rowoffset, [lowcol:highcol] + coloffset) = in;
