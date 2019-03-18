function [ fx fy fz ] = phase_gradient_jc( expphi,dxdydz)
%jclark
%calc gradeint of the phase, without wrapping issues
% see phase_gradient for an alternate method
try
    dxdydz;
catch
    dxdydz=[1,1,1];
end

if numel(dxdydz) == 1,dx=dxdydz;dy=dxdydz;dz=dxdydz;end
if numel(dxdydz) == 2,dx=dxdydz(1);dy=dxdydz(2);dz=1;end
if numel(dxdydz) == 3,dx=dxdydz(1);dy=dxdydz(2);dz=dxdydz(3);end



C=(1/(1i));

if ndims(expphi) == 2
    
    
    [dfx dfy]=gradient(expphi,1);
    
    fx=C*dfx.*conj(expphi);
    fy=C*dfy.*conj(expphi);
end

if ndims(expphi) == 3
   
    [dfx dfy dfz]=gradient(expphi,1);
    fx=C*dfx.*conj(expphi);
    fy=C*dfy.*conj(expphi);
    fz=C*dfz.*conj(expphi);
    
end


end


