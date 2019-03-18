function [ params ] = Matlab_phasing_defaults_1( params,size_data )
%a couple of random defualts
try
    params.pcdi;
catch
    params.pcdi=0;
end

params.sz=size_data;

params.orig=params.pcdi;

try
    params.do=[params.x_not_do,params.y_not_do,params.z_not_do];
catch
    params.do=[0,0,0];
end

try
    params.coh=gauss_3D(11,11,11,params.kernalxy(1),params.kernalxy(2),params.kernalxy(3));
catch
    params.kernalxy=[.1,.1,.1];
    params.coh=gauss_3D(11,11,11,params.kernalxy(1),params.kernalxy(2),params.kernalxy(3));
end
params.update=1;                 %tells the alg weather to update the coh 


end

