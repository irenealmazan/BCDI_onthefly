function [xyz] = cent_arr_autoV2(data,fact)
%centers array by minimising the phase of the autocorrelation, upsamples 
%assumes it is Intensity
%returns the amount required to recenter

try
    fact;
catch
    fact=2;
end

sz=size(data);

newx=fact*sz(2);
newy=fact*sz(1);
if numel(sz) ==3,newz=fact*sz(3);else newz=0;end

%data_up=zero_pad_ver3(data,newx,newy,newz);


%%
xvals=(-50:50)/fact;
yvals=xvals;
zvals=yvals;

auto_o=ifftshift(ifftn(fftshift(data)));

if ndims(data) == 2
    
    [x y]=meshgrid(-fix(sz(2)/2):ceil(sz(2)/2)-1,-fix(sz(1)/2):ceil(sz(1)/2)-1);
    
    S=zero_pad_ver3(ones([size(data)*.25]),size(data,2),size(data,1));
    
    fminx=zeros([1,numel(xvals)]);
    fminy=zeros([1,numel(yvals)]);
    
    for xx=1:numel(xvals);
        
        auto=auto_o.*exp(i*2*pi*xvals(xx).*x/sz(2));
       
        phase=S.*atan2(imag(auto),real(auto));
        
        %fminx(xx)=sum(abs(phase(round(sz(1)/2),:)));
        fminx(xx)=sum(abs(phase(:)));
    end
    
    ind=find(fminx == min(fminx(:)));
    xc=xvals(ind);
    
    for yy=1:numel(yvals);
        
        auto=auto_o.*exp(i*2*pi*yvals(yy).*y/sz(1)).*exp(i*2*pi*xc.*x/sz(2));
       
        phase=S.*atan2(imag(auto),real(auto));
        
        %fminy(yy)=sum(abs(phase(:,round(sz(2)/2))));
        fminy(yy)=sum(abs(phase(:)));
    end        
    
    ind=find(fminy == min(fminy(:)));
    yc=yvals(ind);       
            
    [xyz]=[xc,yc];
    
else
    [x y z]=meshgrid(-fix(sz(2)/2):ceil(sz(2)/2)-1,-fix(sz(1)/2):ceil(sz(1)/2)-1,-fix(sz(3)/2):ceil(sz(3)/2)-1);
    
    S=zero_pad_ver3(ones([size(data)*.25]),size(data,2),size(data,1),size(data,3));
    
    fminx=zeros([1,numel(xvals)]);
    fminy=zeros([1,numel(yvals)]);
    fminz=zeros([1,numel(zvals)]);

    
    for xx=1:numel(xvals);
        
        auto=auto_o.*exp(i*2*pi*xvals(xx).*x/sz(2));
       
        phase=S.*atan2(imag(auto),real(auto));
        
        %fminx(xx)=sum(abs(phase(round(sz(1)/2),:)));
        fminx(xx)=sum(abs(phase(:)));
    end
    
    ind=find(fminx == min(fminx(:)));
    xc=xvals(ind);
    
    for yy=1:numel(yvals);
        
        auto=auto_o.*exp(i*2*pi*yvals(yy).*y/sz(1)).*exp(i*2*pi*xc.*x/sz(2));
       
        phase=S.*atan2(imag(auto),real(auto));
        
        %fminy(yy)=sum(abs(phase(:,round(sz(2)/2))));
        fminy(yy)=sum(abs(phase(:)));
    end        
    
    ind=find(fminy == min(fminy(:)));
    yc=yvals(ind);       
    
    for zz=1:numel(zvals);
        
        auto=auto_o.*exp(i*2*pi*yc.*y/sz(1)).*exp(i*2*pi*xc.*x/sz(2)).*exp(i*2*pi*zvals(zz).*z/sz(3));
       
        phase=S.*atan2(imag(auto),real(auto));
        
        %fminy(yy)=sum(abs(phase(:,round(sz(2)/2))));
        fminz(zz)=sum(abs(phase(:)));
    end        
    
    ind=find(fminz == min(fminz(:)));
    zc=zvals(ind);       
            
    [xyz]=[xc,yc,zc];
   
end
    

end

