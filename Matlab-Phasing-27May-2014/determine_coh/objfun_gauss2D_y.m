function fxy = objfun_gauss2D_y(x,data,modulus,angle)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians



kern = gauss_2D(11,11,0,x(1),angle ) ;
    

guess = convn(modulus.^2,kern,'same');

fxy=sum(sum(  abs(guess-data.^2)    ) )/sum(sum(abs(data.^2)));



end