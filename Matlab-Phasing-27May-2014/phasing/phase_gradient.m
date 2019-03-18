function [ fx fy fz ] = phase_gradient( expphi,dxdydz)
%jclark
%calc gradeint of the phase, without wrapping issues
% see phase_gradient_jc for an alternative method

try
    dxdydz;
catch
    dxdydz=[1,1,1];
end

if numel(dxdydz) == 1,dx=dxdydz;dy=dxdydz;dz=dxdydz;end
if numel(dxdydz) == 2,dx=dxdydz(1);dy=dxdydz(2);dz=1;end
if numel(dxdydz) == 3,dx=dxdydz(1);dy=dxdydz(2);dz=dxdydz(3);end





if ndims(expphi) >= 2
   
    fx=phase_gradient_x(expphi,dx);
    fy=phase_gradient_y(expphi,dy);
    
end

if ndims(expphi) == 3
   
    fz=phase_gradient_z(expphi,dz);
    
end


end

function dphidx = phase_gradient_x(expphi,dx)

dphidx=(1/dx)*angle(circshift(expphi,[0,dx/2,0]).*circshift(conj(expphi),[0,-dx/2,0]));


end

function dphidy = phase_gradient_y(expphi,dy)

dphidy=(1/dy)*angle(circshift(expphi,[dy/2,0,0]).*circshift(conj(expphi),[-dy/2,0,0]));


end


function dphidz = phase_gradient_z(expphi,dz)

dphidz=(1/dz)*angle(circshift(expphi,[0,0,dz/2]).*circshift(conj(expphi),[0,0,-dz/2]));


end