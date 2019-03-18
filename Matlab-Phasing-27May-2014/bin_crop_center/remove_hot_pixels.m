function [data] = remove_hot_pixels(data,thresh)
%jclark

if exist('thresh') ~= 1,thresh=5000;end

switch ndims(data)
    case 4
        for qq=1:size(data,ndims(data))
            for ww=1:size(data,ndims(data)-1)
                L(:,:,ww,qq)=del2(data(:,:,ww,qq));
            end
        end
    case 3
        for qq=1:size(data,ndims(data))
            L(:,:,qq)=del2(data(:,:,qq));
        end
end

ind=find(abs(L) >= thresh);

switch ndims(data)
    case 4
        [aa bb cc dd]=ind2sub(size(L),ind);
        an=aa-1;bn=bb;cn=cc;dn=dd;
    case 3
        [aa bb cc]=ind2sub(size(L),ind);
        an=aa-1;bn=bb;cn=cc;
end
    
%need some indices to use to set the new value to
%choose a neighbouring pixel, make sure it is > 0
if an < 1,an=2;end
if bn < 1,bn=2;end

%now change the values
switch ndims(data)
    case 4
        for qq=1:numel(aa)
            data(aa(qq),bb(qq),cc(qq),dd(qq))=data(an(qq),bn(qq),cn(qq),dn(qq));
        end
    case 3
        for qq=1:numel(aa)
            data(aa(qq),bb(qq),cc(qq))=data(an(qq),bn(qq),cn(qq));
        end
end

end

