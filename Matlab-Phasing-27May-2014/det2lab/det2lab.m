function [pn params] = det2lab(pn,params)
%jclark
%coorindate transformation from detector frame to lab frame
%baed on 34-id-c geometry.  adapted from Ross Harder's 
% python and c code doing the same thing.  See mark pfeifer's thesis
%as well


Nx=size(pn,2);
Ny=size(pn,1);
Nz=size(pn,3);

lambda=params.lam;
delta=params.delta;
gam=params.gam;
dth=params.dth;
dtilt=params.dtilt;
arm=params.arm;
px=params.binning(1)*params.det_px;
py=params.binning(2)*params.det_px;

sx=lambda*arm/px/Nx;
sy=lambda*arm/py/Ny;

if dth ~= 0,sz=lambda/(dth*pi/180)/Nz;else sz=lambda/(dtilt*pi/180)/Nz;end


if isfield(params,'newxyz') == 1
    %resize
    if isempty(params.newxyz) ~= 1
        disp('oooo Resizing sample array.... oooo')
        pn=zero_pad_ver3(pn,params.newxyz(1),params.newxyz(2),params.newxyz(3));
        Nx=params.newxyz(1);Ny=params.newxyz(2);Nz=params.newxyz(3);
        
        %need to change the det size
        px=(arm*lambda/Nx/sx);
        py=(arm*lambda/Ny/sy);
        
        if dth ~= 0,dth=lambda*180/sz/pi/Nz;else dtilt=lambda*180/sz/pi/Nz;end

        
    end
end

use_auto=0;
try
    params.sample_pixel;
    if isempty(params.sample_pixel) == 1,use_auto=1;end
catch
    use_auto=1;
end

if use_auto == 1
    FOVth=0.66;      %maximum reduction allowed in the field of view to prevent cutting the object
    params.sample_pixel=min([sz,sy,sx]);
    params.FOVold=[Nx*sx,Ny*sy,Nz*sz];
    params.FOVnew=[Nx*params.sample_pixel,Ny*params.sample_pixel,Nz*params.sample_pixel];
    FOVratio=min(params.FOVnew./params.FOVold);
    
    if FOVratio(1) < FOVth,params.sample_pixel=FOVth/FOVratio(1)*params.sample_pixel;end
end


if dth ~= 0    
    T=ross_tform(Nx,Ny,Nz,lambda,delta,gam,px,py,dth,arm );
else
    T=ross_tform_dtilt(Nx,Ny,Nz,lambda,delta,gam,px,py,dtilt,arm );
end

params.T=T;
self.size1=Nx;
self.size2=Ny;
self.size3=Nz;

[ self ] = UpdateCoordSystem(self,T(1:3,1:3));

X=reshape(self.coords(1,:),Ny,Nx,Nz);
Y=reshape(self.coords(2,:),Ny,Nx,Nz);
Z=reshape(self.coords(3,:),Ny,Nx,Nz);

X=X-min(X(:));
Y=Y-min(Y(:));
Z=Z-min(Z(:));



disp(['Final sample pixel size (nm) - [',num2str(params.sample_pixel),']'])

XX=self.x*params.sample_pixel*Nx;
YY=self.y*params.sample_pixel*Ny;
ZZ=self.z*params.sample_pixel*Nz;

pn=flipdim(pn,3);

%pre 2012 version
%w = griddata3(X-max(X(:))/2,Y-max(Y(:))/2,Z-max(Z(:))/2,pn,XX-max(XX(:))/2,YY-max(YY(:))/2,ZZ-max(ZZ(:))/2);
%post 2012 version


if verLessThan('matlab', '7.14.0')
    w = griddata3(X-max(X(:))/2,Y-max(Y(:))/2,Z-max(Z(:))/2,pn,XX-max(XX(:))/2,YY-max(YY(:))/2,ZZ-max(ZZ(:))/2);
else
    w = griddata(X-max(X(:))/2,Y-max(Y(:))/2,Z-max(Z(:))/2,pn,XX-max(XX(:))/2,YY-max(YY(:))/2,ZZ-max(ZZ(:))/2);
end



w(isnan(w))=0;

%array=real(w)/max(real(w(:)));
pn=w;

end

function [ self ] = UpdateCoordSystem(self,T)
%%jclark
%adapted from Ross harders c code of the same name
%%
Nx=self.size1;
Ny=self.size2;
Nz=self.size3;
self.dx=1./self.size1;
self.dy=1./self.size2;
self.dz=1./self.size3;

[x y z]=meshgrid( ((1:self.size1)*self.dx),((1:self.size2)*self.dy),(1:self.size3)*self.dz);

self.x=x;
self.y=y;
self.z=z;

r=zeros([3,Ny,Nx,Nz]);

r(1,:,:,:)=x;
r(2,:,:,:)=y;
r(3,:,:,:)=z;

r=reshape(r,3,self.size2*self.size1*self.size3);

%self.coords = nd.dot(r, self.T)
self.coords=T*r;

self.coords=reshape(self.coords,3,self.size2,self.size1,self.size3);
%return self.coords

end

function T=ross_tform(Nx,Ny,Nz,lambda,delta,gam,px,py,dth,arm )
%jclark
%transformation matrix for coordinate transformation
%adapted from Ross Harders c code
%for theta rocking curve

dx=1/Nx;
dy=1/Ny;
dz=1/Nz;

tth=delta;

deg2rad=pi/180;

tth=tth*deg2rad;
gam=gam*deg2rad;
dth=dth*deg2rad;

dpx=px/arm;
dpy=py/arm;

%old pre july 2013
dQdpx(1) = -cos(tth)*cos(gam);
dQdpx(2) = 0.0;
dQdpx(3) = +sin(tth)*cos(gam);

dQdpy(1) = sin(tth)*sin(gam);
dQdpy(2) = -cos(gam);
dQdpy(3) = cos(tth)*sin(gam);

dQdth(1) = -cos(tth)*cos(gam)+1.0;
dQdth(2) = 0.0;
dQdth(3) = sin(tth)*cos(gam);

%new XH version
% dQdpx(1) = cos(tth);
% dQdpx(2) = 0.0;
% dQdpx(3) = -sin(tth);
% 
% dQdpy(1) = sin(tth)*sin(gam);
% dQdpy(2) = -cos(gam);
% dQdpy(3) = cos(tth)*sin(gam);
% 
% dQdth(1) = cos(tth)*cos(gam)-1.0;
% dQdth(2) = 0.0;
% dQdth(3) = -sin(tth)*cos(gam);

%
 
Astar(1) = (2*pi/lambda)*dpx*dQdpx(1);
Astar(2) = (2*pi/lambda)*dpx*dQdpx(2);
Astar(3) = (2*pi/lambda)*dpx*dQdpx(3);

Bstar(1) = (2*pi/lambda)*dpy*dQdpy(1);
Bstar(2) = (2*pi/lambda)*dpy*dQdpy(2);
Bstar(3) = (2*pi/lambda)*dpy*dQdpy(3);

Cstar(1) = (2*pi/lambda)*dth*dQdth(1);
Cstar(2) = (2*pi/lambda)*dth*dQdth(2);
Cstar(3) = (2*pi/lambda)*dth*dQdth(3);

denom=dot(Astar,cross(Bstar,Cstar));
Axdenom=cross(Bstar,Cstar);
Bxdenom=cross(Cstar,Astar);
Cxdenom=cross(Astar,Bstar);

A(1)=2*pi*Axdenom(1)/(denom);
A(2)=2*pi*Axdenom(2)/(denom);
A(3)=2*pi*Axdenom(3)/(denom);

B(1)=2*pi*Bxdenom(1)/(denom);
B(2)=2*pi*Bxdenom(2)/(denom);
B(3)=2*pi*Bxdenom(3)/(denom);

C(1)=2*pi*Cxdenom(1)/(denom);
C(2)=2*pi*Cxdenom(2)/(denom);
C(3)=2*pi*Cxdenom(3)/(denom);


Astarmag=dot(Astar, Astar);
Astarmag = sqrt(Astarmag);
Bstarmag=dot(Bstar, Bstar);
Bstarmag = sqrt(Bstarmag);
Cstarmag=dot(Cstar, Cstar);
Cstarmag = sqrt(Cstarmag);

Amag=dot(A, A);
Amag = sqrt(Amag);
Bmag=dot(B, B);
Bmag = sqrt(Bmag);
Cmag=dot(C, C);
Cmag = sqrt(Cmag);

%sprintf('A* %d  %d  %d',Astar)
%sprintf('B* %d  %d  %d',Bstar)
%sprintf('C* %d  %d  %d',Cstar)

% sprintf('|A*| = %d |A| = %d',[Astarmag,Amag])
% sprintf('|B*| = %d |B| = %d',[Bstarmag,Bmag])
% sprintf('|C*| = %d |C| = %d',[Cstarmag,Cmag])
% 
% sprintf('A %d  %d  %d',A)
% sprintf('B %d  %d  %d',B)
% sprintf('C %d  %d  %d',C)

%sprintf('A.A* %d',dot(A,Astar));
%sprintf('B.B* %d',dot(B,Bstar));
%sprintf('C.C* %d',dot(C,Cstar));

%cpnvert to nanmoeters
%A=A*1e9;
%B=B*1e9;
%C=C*1e9;

T=[A(1) B(1) C(1) 0;A(2) B(2) C(2) 0;A(3) B(3) C(3) 0;0 0 0 1];

end

function T=ross_tform_dtilt(Nx,Ny,Nz,lambda,delta,gam,px,py,dth,arm )
%jclark
%transformation matrix for coordinate transformation
%adapted from Ross Harders c code
%for theta rocking curve

dx=1/Nx;
dy=1/Ny;
dz=1/Nz;

tth=delta;

deg2rad=pi/180;

tth=tth*deg2rad;
gam=gam*deg2rad;
dth=dth*deg2rad;

dpx=px/arm;
dpy=py/arm;


dQdpx(1) = -cos(tth)*cos(gam);
dQdpx(2) = 0.0;
dQdpx(3) = +sin(tth)*cos(gam);

dQdpy(1) = sin(tth)*sin(gam);
dQdpy(2) = -cos(gam);
dQdpy(3) = cos(tth)*sin(gam);

dQdth(1) = 0.0;
dQdth(2) = cos(tth)*cos(gam)-1.0;
dQdth(3) = -sin(gam);

 
Astar(1) = (2*pi/lambda)*dpx*dQdpx(1);
Astar(2) = (2*pi/lambda)*dpx*dQdpx(2);
Astar(3) = (2*pi/lambda)*dpx*dQdpx(3);

Bstar(1) = (2*pi/lambda)*dpy*dQdpy(1);
Bstar(2) = (2*pi/lambda)*dpy*dQdpy(2);
Bstar(3) = (2*pi/lambda)*dpy*dQdpy(3);

Cstar(1) = (2*pi/lambda)*dth*dQdth(1);
Cstar(2) = (2*pi/lambda)*dth*dQdth(2);
Cstar(3) = (2*pi/lambda)*dth*dQdth(3);

denom=dot(Astar,cross(Bstar,Cstar));
Axdenom=cross(Bstar,Cstar);
Bxdenom=cross(Cstar,Astar);
Cxdenom=cross(Astar,Bstar);

A(1)=2*pi*Axdenom(1)/(denom);
A(2)=2*pi*Axdenom(2)/(denom);
A(3)=2*pi*Axdenom(3)/(denom);

B(1)=2*pi*Bxdenom(1)/(denom);
B(2)=2*pi*Bxdenom(2)/(denom);
B(3)=2*pi*Bxdenom(3)/(denom);

C(1)=2*pi*Cxdenom(1)/(denom);
C(2)=2*pi*Cxdenom(2)/(denom);
C(3)=2*pi*Cxdenom(3)/(denom);


Astarmag=dot(Astar, Astar);
Astarmag = sqrt(Astarmag);
Bstarmag=dot(Bstar, Bstar);
Bstarmag = sqrt(Bstarmag);
Cstarmag=dot(Cstar, Cstar);
Cstarmag = sqrt(Cstarmag);

Amag=dot(A, A);
Amag = sqrt(Amag);
Bmag=dot(B, B);
Bmag = sqrt(Bmag);
Cmag=dot(C, C);
Cmag = sqrt(Cmag);

%sprintf('A* %d  %d  %d',Astar)
%sprintf('B* %d  %d  %d',Bstar)
%sprintf('C* %d  %d  %d',Cstar)

% sprintf('|A*| = %d |A| = %d',[Astarmag,Amag])
% sprintf('|B*| = %d |B| = %d',[Bstarmag,Bmag])
% sprintf('|C*| = %d |C| = %d',[Cstarmag,Cmag])
% 
% sprintf('A %d  %d  %d',A)
% sprintf('B %d  %d  %d',B)
% sprintf('C %d  %d  %d',C)

%sprintf('A.A* %d',dot(A,Astar));
%sprintf('B.B* %d',dot(B,Bstar));
%sprintf('C.C* %d',dot(C,Cstar));

%cpnvert to nanmoeters
%A=A*1e9;
%B=B*1e9;
%C=C*1e9;

T=[A(1) B(1) C(1) 0;A(2) B(2) C(2) 0;A(3) B(3) C(3) 0;0 0 0 1];

end