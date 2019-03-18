function fxy = objfun_standarddev(x,array1,array2,array3)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


xtot=sum(x(:));
x=x/xtot;           %need it to be normlaized

x1=x(1);
x2=x(2);
x3=x(3);

arr=x1*abs(array1)+x2*abs(array2)+x3*abs(array3);


fxy=std(arr(:).^4)*1e12;
arr=0;

end
