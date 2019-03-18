function [ threshold ] = calc_percent_from_auto( data)
%jclark

auto_c=abs(fftshift(ifftn(fftshift(data.^2))));
temp=shrink_wrap(abs(auto_c).^.5,.1,.2);  %.^.5,.05,.1
threshold=sum(temp(:))/prod(size(data));

if ndims(data) == 3,threshold=threshold/8.0;end
if ndims(data) == 2,threshold=threshold/4.0;end

clear temp
clear auto_c


end

