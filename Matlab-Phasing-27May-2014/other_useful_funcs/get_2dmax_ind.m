function [x y I0] = get_2dmax_ind(data)
%jclark

[I0 v]=max(data(:));
[y x]=ind2sub(size(data),v);


end

