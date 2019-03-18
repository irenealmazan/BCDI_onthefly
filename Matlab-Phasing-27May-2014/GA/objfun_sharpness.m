function fxy = objfun_sharpness(x,array1,array2,array3,norm)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


xtot=sum(x(:));
x=x/xtot;           %need it to be normlaized

x1=x(1);
x2=x(2);
x3=x(3);

arr=x1*abs(array1)+x2*abs(array2)+x3*abs(array3);
%arr=x1*(array1)+x2*(array2)+x3*(array3);


fxy=sum(abs(arr(:)).^4)/norm;
arr=0;

end

