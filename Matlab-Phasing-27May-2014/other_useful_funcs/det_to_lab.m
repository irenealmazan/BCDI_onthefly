function [ Tarray ] = det_to_lab(pn,delta,gamma,dpx,dpy,dt,arm,lambda )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

sz=size(pn);
Nx=sz(2);
Ny=sz(1);
Nt=sz(3);

%pn=rot3D(pn,1);
%sz=size(pn);
%for qq=1:sz(3),pn(:,:,qq)=flipud(pn(:,:,qq));end
%pn=flipdim(pn,3);

%Tr=transform_matrix(delta,gamma,Nx,Ny,Nt,dpx,dpy,dt,arm,lambda);
Tr=rotation_matrix3d(delta,gamma,0);
%%Tr=ross_tform(Nx,Ny,Nt,lambda,delta,gamma,dpx,dpy,dt,arm);
Tr=Tr/max(abs(Tr(:)));
Tr(4,4)=1;

center=(size(pn)+1)/2;

T1=[1 0 0 0;0 1 0 0;0 0 1 0;-center 1;];

%Tr=[cos(delta) 0 -sin(delta) 0;0 1 0 0;sin(delta) 0 cos(delta) 0;0 0 0 1];

T2=[1 0 0 0;0 1 0 0;0 0 1 0;center 1;];

T=T1*Tr*T2;

tform=maketform('affine',T);

R=makeresampler('linear','fill');

tdims_a=[1 2 3];
tdims_b=[1 2 3];
tsize_b=[Ny Nx Nt];
tmap_b=[];
F=0;

Tarray=tformarray(pn,tform,R,tdims_a,tdims_b,tsize_b,tmap_b,F);


end

