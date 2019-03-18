function [ support ] = support_from_phase(phase)
%jclark
%get support from objects phase based on gradient and  smoothing operations 

I=phase;

[~, threshold] = edge(I, 'sobel');
fudgeFactor = .25;
BWs = edge(I,'sobel', threshold * fudgeFactor);

%se90 = strel('line', 3, 90);
%se0 = strel('line', 3, 0);
%BWsdil = imdilate((BWs), [se90 se0]);

ngrowth=20;
BWsdil=BWs;
for qq=1:ngrowth
    BWsdil=BWsdil+grow_support(BWs,1.0+qq*.005);
end
BWsdil=abs(1-ceil(BWsdil/max(BWsdil(:))));


BWdfill = imfill(BWsdil, 'holes');

seD = strel('diamond',1);
BWfinal = imerode(BWdfill,seD);
BWfinal = imerode(BWfinal,seD);

support=BWfinal;

end

