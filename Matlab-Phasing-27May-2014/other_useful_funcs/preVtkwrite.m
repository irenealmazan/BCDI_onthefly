function [ NewCoords ] = preVtkwrite(Nx, Ny, Nz,arm, lam, delta, gam, dpx, dpy, dth,dx ,dy,dz)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%dx,dy,dz is 1/Nx etc
%dpx = px/arm
size1=Nx;
size2=Ny;
size3=Nz;

[ T ] = transform_matrix(delta,gam,Nx,Ny,Nz,dpx,dpy,dth,arm,lam );

for ix=1:size1
  for iy=1:size2
    for iz=1:size3
      
        index = 3*(iz+size3*(iy+size2*ix));

        Rold(1)=(size1-ix-1)*dx;
        Rold(2)=iy*dy;
        Rold(3)=iz*dz;
        %mytransform(Rold, Rnew, T);
        Rnew = Rot_tform( Rold,T );
        data(index+0) = Rnew(1);
        data(index+1) = Rnew(2);
        data(index+2) = Rnew(3);
    end
  end
end



end

