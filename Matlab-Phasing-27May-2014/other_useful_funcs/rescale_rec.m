function [ array ] = rescale_rec(pn,dpx,dpy,dt,arm,lambda )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
deg2rad=pi/180.0;

[Ny Nx Nz]=size(pn);
N0=[Nx, Ny, Nz];
Nx0=N0(1);
Ny0=N0(2);
Nz0=N0(3);

dpt=sin(dt*deg2rad*Nz0)/Nz0;%arm*tan(dt*deg2rad);

dsx0=arm*lambda/Nx/dpx;
dsy0=arm*lambda/Ny/dpy;
dsz0=arm*lambda/Nz/dpt;

%xang=atan(Nx/2*dpx/arm)*2;
%yang=atan(Ny/2*dpy/arm)*2;

pxs=[dsx0,dsy0,dsz0];

minpx=find(pxs == min(pxs));

tol=0.05;

switch minpx
    
    case 1;
        Ny=Nx*dpx/dpy;
        Nz=Nx*dpx/dpt;
    case 2;
        Nx=Ny*dpy/dpx;
        Nz=Ny*dpy/dpt;
    case 3;
        Nx=Nz*dpt/dpx;
        Ny=Nz*dpt/dpy;
end

N1=round([Nx Ny Nz]);
Nx=N1(1);
Ny=N1(2);
Nz=N1(3);

dsx1=arm*lambda/Nx/dpx;
dsy1=arm*lambda/Ny/dpy;
dsz1=arm*lambda/Nz/dpt;

nnc=[0,0,0,0,0,0];

fftpn=ifftshift(fftn(fftshift(pn)));

%if mod(Nx,2) == 0, nnc(1:2)=abs(Nx0-Nx)/2;else
    %nnc(1:2)=[floor
nnc=[floor(abs(Nx0-Nx)/2),ceil(abs(Nx0-Nx)/2),floor(abs(Ny0-Ny)/2), ...
    ceil(abs(Ny0-Ny)/2),floor(abs(Nz0-Nz)/2),ceil(abs(Nz0-Nz)/2)];

fftpn=init_pad(fftpn,nnc);
pn0=ifftshift(ifftn(fftshift(fftpn)));
pn0=center_array(pn0);
pn=pn0;
%pn=zero_phase(pn0);
%pn=remove_ramp(pn);

xyz=center_of_mass(abs(pn));
xyz=-1*round([xyz(2),xyz(1),xyz(3)]);
pn=circshift(real(pn),xyz)+1i*circshift(imag(pn),xyz);

NN=max(size(pn));
%pn=zero_pad_ver2(pn,N0(1),N0(2),N0(3));
%pn=zero_pad(pn,N0(1),N0(2),N0(3));
%pn=zero_pad_ver2(pn,NN,NN,NN);
array=pn;
%delta=s1;    
end

function data = init_pad(data,nnc)

if numel(nnc) == 6
    
    if sum(nnc(1:2)) > 0, 
        disp('doing intial x padding....')
        data=padarray(data,[0 abs(nnc(1)) 0],0,'pre');
        data=padarray(data,[0 abs(nnc(2)) 0],0,'post');
    end
    
    if sum(nnc(3:4)) > 0, 
        disp('doing intial y padding....')
        data=padarray(data,[abs(nnc(3)) 0 0],0,'pre');
        data=padarray(data,[abs(nnc(4)) 0 0],0,'post');
    end
    
    if sum(nnc(5:6)) > 0, 
        disp('doing intial z padding....')
        data=padarray(data,[0 0 abs(nnc(5))],0,'pre');
        data=padarray(data,[0 0 abs(nnc(6))],0,'post');
    end
    
else
    disp('if initial pading is required')
    disp('set nnc=[+x0,+x1,+y0,+y1,+z0,+z1], where nnc gives pixels to pad')
    disp('onto each dimension otherwise set nnc=[0]')
end


end

function pn = zero_phase(pn)

%set COM phase to zero, use as a reference
disp('removing phase offset....')
disp('')
sz=size(pn);
i=round(sz(1)/2);
j=round(sz(2)/2);
k=round(sz(3)/2);
phi0=atan2(imag(pn(i,j,k)),real(pn(i,j,k)));
pn=pn*exp(-1i*phi0);

end

function pn = remove_ramp(pn)

disp('removing phase ramp....')
disp('')
amp=abs(pn);
sz=size(pn);
i=round(sz(1)/2);
j=round(sz(2)/2);
k=round(sz(3)/2);
damp=fftshift(fftn(fftshift(amp)));
[l,m,n]=ind2sub(size(damp),find( abs(damp) == max(max(max(abs(damp))))));
dpn=circshift(fftshift(fftn(fftshift(pn))),[l-i,m-j,n-k]);
pn=fftshift(ifftn(fftshift(dpn)));

end