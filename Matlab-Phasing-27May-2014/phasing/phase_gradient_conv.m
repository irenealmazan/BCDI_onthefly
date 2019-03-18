function [fx hx fy hy] = phase_gradient_conv(expphi)
%jclark
%calc gradeint of the phase, without wrapping issues
% see phase_gradient_jc for an alternative method

C=(1/(1i));

if ndims(expphi) >= 2
   
    [dfx hx]=phase_gradient_x(expphi);
    [dfy hy]=phase_gradient_y(expphi);
    fx=C*dfx.*conj(expphi);
    fy=C*dfy.*conj(expphi);
end

if ndims(expphi) == 3
   
    [dfz hz]=phase_gradient_z(expphi);
    fx=C*dfx.*conj(expphi);
    fy=C*dfy.*conj(expphi);
    fz=C*dfz.*conj(expphi);
end


end

function [dphidx h]= phase_gradient_x(expphi,dxt)

h=fspecial('sobel');
sz=size(expphi);
h=h/sum(abs(h(:)));
if ndims(expphi) == 2
    h=fftshift(fftn(fftshift(zero_pad_ver3(h',sz(2),sz(1)))));
end
if ndims(expphi) == 3
    h=fftshift(fftn(fftshift(zero_pad_ver3(h',sz(2),sz(1),sz(3)))));
end

dphidx=fftshift(ifftn(fftshift(fftshift(fftn(fftshift(expphi))).*h)));


end

function [dphidy h] = phase_gradient_y(expphi)

h=fspecial('sobel');
sz=size(expphi);
h=h/sum(abs(h(:)));
if ndims(expphi) == 2
    h=fftshift(fftn(fftshift(zero_pad_ver3(h,sz(2),sz(1)))));
    
end
if ndims(expphi) == 3
    h=fftshift(fftn(fftshift(zero_pad_ver3(h,sz(2),sz(1),sz(3)))));
    
end

dphidy=fftshift(ifftn(fftshift(fftshift(fftn(fftshift(expphi))).*h)));



end


function dphidz = phase_gradient_z(expphi)

h=fspecial('sobel');
sz=size(expphi);

if ndims(expphi) == 2
    h=zero_pad_ver3(h',sz(2),sz(1));
end
if ndims(expphi) == 3
    h=zero_pad_ver3(h',sz(2),sz(1),sz(3));
end

dphidy=fftshift(fftn(fftshift(expphi.*h)));



end