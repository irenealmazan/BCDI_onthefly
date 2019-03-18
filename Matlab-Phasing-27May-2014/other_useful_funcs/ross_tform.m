function T=ross_tform(Nx,Ny,Nz,lambda,delta,gam,px,py,dth,arm )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
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

dQdth(1) = -cos(tth)*cos(gam)+1.0;
dQdth(2) = 0.0;
dQdth(3) = sin(tth)*cos(gam);

  

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

sprintf('|A*| = %d |A| = %d',[Astarmag,Amag])
sprintf('|B*| = %d |B| = %d',[Bstarmag,Bmag])
sprintf('|C*| = %d |C| = %d',[Cstarmag,Cmag])

sprintf('A %d  %d  %d',A)
sprintf('B %d  %d  %d',B)
sprintf('C %d  %d  %d',C)

%sprintf('A.A* %d',dot(A,Astar));
%sprintf('B.B* %d',dot(B,Bstar));
%sprintf('C.C* %d',dot(C,Cstar));

%cpnvert to nanmoeters
%A=A*1e9;
%B=B*1e9;
%C=C*1e9;


T=[A(1) B(1) C(1) 0;A(2) B(2) C(2) 0;A(3) B(3) C(3) 0;0 0 0 1];

end

