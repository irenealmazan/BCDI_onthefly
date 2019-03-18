function [xyz] = cent_arr_auto(data,fact)
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

data_up=zero_pad_ver3(data,newx,newy,newz);

%%
xvals=-5:5;
yvals=xvals;
zvals=yvals;



if ndims(data) == 2
    
    S=zero_pad_ver3(ones([size(data_up)*.25]),size(data_up,2),size(data_up,1));
    
    fminx=zeros([1,numel(xvals)]);
    fminy=zeros([1,numel(yvals)]);
    
    for xx=1:numel(xvals);
        auto=fftshift(ifftn(fftshift(circshift(data_up,[0,xvals(xx)]))));
        %if xx == 1,S=(abs(auto) > .05 * max( abs(auto(:))));end
        phase=S.*atan2(imag(auto),real(auto));
        fminx(xx)=sum(abs(phase(:)));
    end
    ind=find(fminx == min(fminx(:)));
    xc=xvals(ind);
    for yy=1:numel(yvals);
        auto=fftshift(ifftn(fftshift(circshift(data_up,[yvals(yy),xc])))); 
        %if yy == 1,S=(abs(auto) > .05 * max( abs(auto(:))));end
        phase=S.*atan2(imag(auto),real(auto));
        fminy(yy)=sum(abs(phase(:)));
    end        
    
    ind=find(fminy == min(fminy(:)));
    yc=yvals(ind);       
            
    [xyz]=[xc,yc]/fact;
    
else
    S=zero_pad_ver3(ones([size(data_up)*.25]),size(data_up,2),size(data_up,1),size(data_up,3));
    fminx=zeros([1,numel(xvals)]);
    fminy=zeros([1,numel(yvals)]);
    fminz=zeros([1,numel(zvals)]);    
    
    for xx=1:numel(xvals);
        auto=fftshift(ifftn(fftshift(circshift(data_up,[0,xvals(xx),0]))));
        %if xx == 1,S=(abs(auto) > .05 * max( abs(auto(:))));end
        phase=S.*atan2(imag(auto),real(auto));
        %fminx(xx)=sum(abs(phase(:,:,:)));
        fminx(xx)=sum(abs(phase(round(sz(1)/2),:,round(sz(3)/2))));
    end
    ind=find(fminx == min(fminx(:)));
    xc=xvals(ind(1));
    for yy=1:numel(yvals);
        auto=fftshift(ifftn(fftshift(circshift(data_up,[yvals(yy),xc,0])))); 
        %if yy == 1,S=(abs(auto) > .05 * max( abs(auto(:))));end
        phase=S.*atan2(imag(auto),real(auto));
        %fminy(yy)=sum(abs(phase(:)));
        fminy(yy)=sum(abs(phase(:,round(sz(2)/2),round(sz(3)/2))));
    end        
    
    ind=find(fminy == min(fminy(:)));
    yc=yvals(ind(1));       
         
    
    for zz=1:numel(yvals);
        auto=fftshift(ifftn(fftshift(circshift(data_up,[yc,xc,zvals(zz)]))));
        %if zz == 1,S=(abs(auto) > .05 * max( abs(auto(:))));end
        phase=S.*atan2(imag(auto),real(auto));
        %fminz(zz)=sum(abs(phase(:)));
        fminz(zz)=sum(abs(phase(round(sz(1)/2),round(sz(2)/2),:)));
    end
    ind=find(fminz == min(fminz(:)));
    zc=zvals(ind(1));
    
    [xyz]=[xc,yc,zc]/fact;
   
end
    
    
%%
% lb=[-5,-5,-5];
% ub=[5,5,5];
% 
% options = optimset('Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-6);
% %options=psoptimset('Display','off','TolFun',1e-7,'TolMesh',1e-6,'MeshContraction',.5,'PollingOrder','Random','MeshExpansion',3);
%    
% if numel(sz) == 3,x0=[.5,.5,.5];else x0=[.5,.5];end
% 
% f=@(x)objfun_min_phase(x,data);
% 
% %x=patternsearch(f,x0,[],[],[],[],[lb],[ub],[],options);
% [x fval]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
% 
% disp([x])



end

function fxy = objfun_min_phase(x,array)
%array is the intensity
%x=[x,y,z]
x=round(x);

row_shift=x(1);
col_shift=x(2);

if numel(x) == 3,z_shift=x(3);end

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% buf2ft=fftn(array);
% [nr,nc,nz]=size(buf2ft);
% Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
% Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
% 
% if ndims(array) == 3
%     Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
%     [Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
%     Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
% else
%     [Nc,Nr] = meshgrid(Nc,Nr);
%     Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
% end
%%
if numel(x) == 3
    Greg=ifftshift(ifftn(fftshift(circshift(array,[col_shift, row_shift, z_shift]))));
else
    Greg=ifftshift(ifftn(fftshift(circshift(array,[col_shift, row_shift]))));
end

%%
sz=size(Greg);

if numel(sz) == 3    
    roi=Greg( sz(1)/2-sz(1)/4:sz(1)/2+sz(1)/4,sz(2)/2-sz(2)/4:sz(2)/2+sz(2)/4,sz(3)/2-sz(3)/4:sz(3)/2+sz(3)/4);       
else
    roi=Greg( sz(1)/2-sz(1)/4:sz(1)/2+sz(1)/4,sz(2)/2-sz(2)/4:sz(2)/2+sz(2)/4);   
end

fxy=sum( abs( atan2(imag(roi(:)), real(roi(:)))));

end

% if ndims(data) == 2
%    
%     fmin=zeros([numel(yvals),numel(xvals)]);
%     xarr=zeros([numel(yvals),numel(xvals)]);
%     yarr=zeros([numel(yvals),numel(xvals)]);
%     
%     for xx=1:numel(xvals);
%         for yy=1:numel(yvals);
%             
%             
%             auto=fftshift(ifftn(fftshift(circshift(data,[yvals(yy),xvals(xx)]))));
%             
%             if sum(xx-1+yy-1) ==0,S=(abs(auto) > .05 * max( abs(auto(:))));end
%             
%             phase=S.*atan2(imag(auto),real(auto));
%             fmin(yy,xx)=sum(abs(phase(:)));
%             
%             xarr(yy,xx)=xvals(xx);
%             yarr(yy,xx)=yvals(yy);
%             
%         end
%     end
%     
%     ind=( fmin == min(fmin(:)));
%     
%     [xyz]=[xarr(ind),yarr(ind)]/fact;
%     
% else
%     
%     fmin=zeros([numel(yvals),numel(xvals),numel(zvals)]);
%     xarr=zeros([numel(yvals),numel(xvals),numel(zvals)]);
%     yarr=zeros([numel(yvals),numel(xvals),numel(zvals)]);
%     zarr=zeros([numel(yvals),numel(xvals),numel(zvals)]);
%     
%     
%     for xx=1:numel(xvals);
%         for yy=1:numel(yvals);
%             for zz=1:numel(zvals);
% 
%                 auto=fftshift(ifftn(fftshift(circshift(data_up,[yvals(yy),xvals(xx),zvals(zz)]))));
% 
%                 if sum(xx-1+yy-1+zz-1) ==0,S=(abs(auto) > .05 * max( abs(auto(:))));end
% 
%                 phase=S.*atan2(imag(auto),real(auto));
%                 fmin(yy,xx,zz)=sum(abs(phase(:)));
%                 xarr(yy,xx)=xvals(xx);
%                 yarr(yy,xx)=yvals(yy);
%                 zarr(yy,xx)=zvals(zz);
%                 
%             end
%         end
%     end
%     
%     ind=( fmin == min(fmin(:)));
%     [xyz]=[xarr(ind),yarr(ind),zarr(ind)]/fact;
% end