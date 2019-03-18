function [ circ ] = generate_circle_nd(nn,rad,dim)
%jclark

if exist('dim') ~= 1 ,dim=2;end

if dim == 2
    [x y]=meshgrid(-nn(1)/2+1:nn(1)/2,-nn(2)/2+1:nn(2)/2 );
    z=0;
end

if dim == 3
    [x y z]=meshgrid(-nn(1)/2+1:nn(1)/2,-nn(2)/2+1:nn(2)/2,-nn(3)/2+1:nn(3)/2);
end


circ=( (x.^2+y.^2+z.^2).^.5 <=rad )*1e0;




end

