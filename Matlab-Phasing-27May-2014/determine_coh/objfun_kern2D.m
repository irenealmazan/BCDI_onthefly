function fxy = objfun_kern2D(x,data,modulus)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians


kern = x;

kern=kern/sum(sum(sum(kern)));

guess = convn(modulus.^2,kern,'same');

fxy=sum(sum(sum(  abs(guess-data.^2)    ) )) /sum(sum(sum(abs(data.^2))));



end