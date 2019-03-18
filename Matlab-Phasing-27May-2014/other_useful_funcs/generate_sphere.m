function [ circ ] = generate_sphere(nn,rad )
%jclark

if mod(nn,2) == 1
    [X Y Z]=meshgrid(-nn/2+0.5:nn/2-0.5);
end

if mod(nn,2) == 0
    [X Y Z]=meshgrid( -(nn-1)/2:(nn-1)/2 );
end


circ=( (X.^2+Y.^2+Z.^2).^.5 <= rad )*1e0;




end

