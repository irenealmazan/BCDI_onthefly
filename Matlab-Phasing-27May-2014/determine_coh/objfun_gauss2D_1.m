function fxy = objfun_gauss2D_1(x,data,modulus,do)
%2d gauss, checks for even or odd length so as to center it correctly and
%avoid division by a 0 sigma.  0 sigma can be entered and it will return
%the correct gauss form, i.e repeated 1d gaussians

ind=find(do(1:2) == 0); %get the one to do

if ind == 1,
    x(2)=0;
    
    if do(3) == 0
        angle=x(3);
        kern = gauss_2D(11,11,x(1),0,angle ) ;
    else 
        x(3)=0;
        angle=x(3);
        kern = gauss_2D(11,11,x(1),0,angle ) ;
end

if ind == 2,
    x(1)=0;
    
    if do(3) ==0,
        angle=x(3);
        kern = gauss_2D(11,11,0,x(2),angle ) ;
    else
        x(3)=0;
        angle=x(3);
        kern = gauss_2D(11,11,0,x(2),angle);
    end
end

guess = convn(modulus.^2,kern,'same');

fxy=sum(sum(  abs(guess-data.^2)    ) )/sum(sum(abs(data.^2)));



end