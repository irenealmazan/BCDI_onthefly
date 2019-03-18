function [ trigger ] = alg_trigger(int1,int2,iterations)
%jclark

trigger=sort([(0:(int1+int2):iterations),((0:(int1+int2):iterations)+int1)]);


end

