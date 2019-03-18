function [ image ] = grey_level_my_image( image, levels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

amp=abs(image);
tot=sum(amp(:));

ph=angle(image);

amp=round(amp/max(amp(:))*levels);

amp=amp/sum(amp(:))*tot;

image=amp.*exp(i*ph);

end

