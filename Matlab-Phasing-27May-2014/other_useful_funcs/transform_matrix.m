function [ T ] = transform_matrix(delta,gamma,Nx,Ny,Nt,dpx,dpy,dt,arm,lambda )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
deg2rad=pi/180e0;

delta=delta*deg2rad;
gamma=gamma*deg2rad;

dt=dt*deg2rad;
dpx=dpx/arm;
dpy=dpy/arm;

AA=1e0/( sin(delta)*cos(gamma))*lambda/2.0/pi;

Sx=1/Nx/dpx*AA;
Sy=1/Ny/dpy*AA;
St=1/Nt/dt*AA;

A1=-sin(delta)*cos(gamma)^2*Sx;
A2=( cos(delta)*sin(gamma)-sin(gamma)*cos(gamma) )*Sx;
A3=(cos(gamma)-cos(delta)*cos(gamma)^2)*Sx;

B1=0;
B2=-sin(delta)*Sy;
B3=0;

C1=sin(delta)*cos(gamma)*St;
C2=sin(gamma)*St;
C3=cos(delta)*cos(gamma)*St;

T=[A1 B1 C1 0;A2 B2 C2 0;A3 B3 C3 0;0 0 0 1];
%T=T/max(max(abs(T)));
%T(4,4)=1;
% %T=[0.9239         0   -0.3827         0 ...
%         ; 0    1.0000         0         0 ...
%    ; 0.3827         0    0.9239         0  ...
%   ;-12.5691         0   18.8110    1.0000  ]



end

