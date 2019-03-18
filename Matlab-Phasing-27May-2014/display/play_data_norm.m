function  play_data_norm( data,orient,cplx )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
try
    orient;
    if isempty(orient),orient='xy';end
catch
    orient='xy';
end

try
    cplx;
catch
    cplx=0;
    if isreal(data) == 0,cplx=1;end
end

figure(73)

switch orient
        
        case {'xy','yx'}
            nn=size(data,3);
        case {'zx','xz'}
            nn=size(data,1);
        case {'zy','yz'}
            nn=size(data,2);
    
end

for qq=1:nn
    title(num2str(qq))
    pause(.1)

    if cplx == 0
        switch orient

            case {'xy','yx'}
                imagesc((data(:,:,qq))); 
                %imagesc((data(:,:,qq)))
            case {'zx','xz'}
                imagesc(squeeze(data(qq,:,:)))
            case {'zy','yz'}
              imagesc(squeeze(data(:,qq,:)))

        end
    else
        switch orient
        
            case {'xy','yx'}
                imagesc(c2image(data(:,:,qq)))
            case {'zx','xz'}
                imagesc(c2image(squeeze(data(qq,:,:))))
            case {'zy','yz'}
              imagesc(c2image(squeeze(data(:,qq,:))))

        end
 
    end
end


end

