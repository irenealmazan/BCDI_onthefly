function [U S V isv rsv] = svd_my_images( array)
%jclark
%assumes time is the last dim
sz=size(array);

array=reshape(array,[prod(sz(1:end-1)),sz(end)]);  %reshape the vectors

[U S V]=svd(array,0);

isv=reshape(U*S,[sz]);
rsv=S*V;

end

