function [ TV ] = calc_TV(array)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if ndims(array) == 2
    
    [Fx Fy]=gradient(array);
    TV=abs(Fx)+abs(Fy);
end


if ndims(array) == 3
    
    [Fx Fy Fz]=gradient(array);

    TV=abs(Fx)+abs(Fy)+abs(Fz);
end


end

