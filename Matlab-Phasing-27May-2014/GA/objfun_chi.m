function fxy = objfun_chi(x,array1,array2,array3,data)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


xtot=sum(x(:));
x=x/xtot;           %need it to be normlaized

x1=x(1);
x2=x(2);
x3=x(3);

ph_alpha=atan2(imag(array1),real(array1));
arr=x1*array1+x2*array2;%(x1*abs(array1)+x2*abs(array2)+x3*abs(array3)).*exp(i*ph_alpha);

arr=fftshift(fftn(fftshift(arr)));

fxy=calc_chi(abs(arr),data);

arr=0;

end