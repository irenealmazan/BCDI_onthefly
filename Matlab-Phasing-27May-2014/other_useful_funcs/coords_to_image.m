function [ image] = coords_to_image(x,y,nx,ny)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

image=zeros([ny nx]);

for qq=1:numel(x)

    
    xc=round(x(qq))+round(nx/2);
    yc=round(y(qq))+round(ny/2);
    try
        image(yc,xc)=image(yc,xc)+1;
    catch
    end
end


end

