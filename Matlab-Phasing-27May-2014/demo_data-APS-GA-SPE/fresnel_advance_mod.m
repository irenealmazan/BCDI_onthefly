function [eq12] = fresnel_advance_mod(U0, dx, dy, z, lambda)
% The function receives a field U0 at wavelength lambda
% and returns the field U after distance z, using the Fresnel
% approximation. dx, dy, are spatial resolution.

%the equation numbers refer to the equation numbers in this paper:
%Mas, D., Garcia, J., Ferreira, C., Bernardo, L. M. & Marinho, F. Fast algorithms for free-space diffraction patterns calculation. Opt. Commun. 164, 233?245 (1999).

k=2*pi/lambda;
[ny, nx] = size(U0); 

Lx = dx * nx;
Ly = dy * ny;
% 
dfx = 1./Lx;
dfy = 1./Ly;

% u = ones(nx,1)*((1:nx)-nx/2)*dfx;    
% v = ((1:ny)-ny/2)'*ones(1,ny)*dfy;   
% 
% O = fftshift(fft2(U0));
% 
% H = exp(1i*k*z)*exp(-1i*pi*lambda*z*(u.^2+v.^2));  %original expression
% 
% U = ifft2(ifftshift(O.*H));  

sz = size(U0);

[X,Y ] = meshgrid([-sz(2)/2:sz(2)/2-1], [-sz(1)/2:sz(1)/2-1]);

%dx is the pixel size
%  eq6 = fftshift(exp(1i*pi*lambda*z*(X.^2+Y.^2)/(sz(1)*dx)^2)).*...
%      fftshift(fft2(U0.*exp(1i*pi*(sz(1)*dx)^2*(X.^2 + Y.^2)/(lambda*z*sz(1)^2))));


H2 = exp(-1i*pi*lambda*z*(X.^2+Y.^2)/(sz(1)*dx)^2);
O2 = fftshift(fft2(U0));
eq12 = ifft2(ifftshift(O2.*H2)); 


% figure(6); imagesc(angle(eq6)); title('eq6');
% figure(7); imagesc(angle(U)); title('original');
% figure(8); imagesc(angle(eq12)); title('eq12');
% 
% figure(1); imagesc(abs(eq6)); title('eq6');
% figure(2); imagesc(abs(U)); title('original');
% figure(3); imagesc(abs(eq12)); title('eq12');
% 
% figure(1); imagesc(fftshift(abs(fft2(eq6)))); title('eq6');
% figure(2); imagesc(fftshift(abs(fft2(U)))); title('original');
% figure(3); imagesc(fftshift(abs(fft2(eq12)))); title('eq12');

