function fxy = least_squares(x,data,modulus)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


kern = abs(x) ;

kern=kern/sum(sum(sum(kern)));

guess = convn(modulus.^2,kern,'same');

fxy=abs(guess-data.^2)/sum(sum(sum(abs(data.^2))));

end

