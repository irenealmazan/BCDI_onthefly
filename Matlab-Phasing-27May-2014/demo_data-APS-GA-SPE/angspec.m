%angular spectrum method propogation example

clc; clear; close all; % parameters 
layer = zeros(499, 499); diameter = 100e-6; %meters
radius = diameter/2; lambda = 500e-9; % wavelength (meters) 
k = 2*pi/lambda; % wavenumber 
z = 2e-4; % distance (meters) 
phy_x = 1e-3; % physical width (meters) 
phy_y = 1e-3; % physical length (meters) 
obj_size = size(layer); 

% generate meshgrid 
dx = linspace(-phy_x/2, phy_x/2, obj_size(2)); dx = dx(1:length(dx)); 
dy = linspace(-phy_y/2, phy_y/2, obj_size(1)); dy = dy(1:length(dy)); 

% generate a circle aperture 
for i = 1:obj_size(2) 
    for j = 1:obj_size(1) 
        if sqrt(dx(i)^2 + dy(j)^2) <= radius; layer(j, i) = 1; 
        end
    end
end

Fs_x = obj_size(2)/phy_x; Fs_y = obj_size(1)/phy_y; dx2 = Fs_x^(-1); dy2 = Fs_y^(-1); 
dFx = Fs_x/obj_size(2); dFy = Fs_y/obj_size(1); 
Fx = (-Fs_x/2:dFx:(Fs_x/2 - dFx)); Fy = (-Fs_y/2:dFy:(Fs_y/2 - dFy)); 
% alpha and beta (wavenumber components) 
alpha = lambda.*Fx; beta = lambda.*Fy; 
% gamma_cust 
gamma_cust = zeros(length(beta), length(alpha)); 
for j = 1:length(beta) 
    for i = 1:length(alpha) 
        if (alpha(i)^2 + beta(j)^2) > 1 gamma_cust(j, i) = 0; 
        else gamma_cust(j, i) = sqrt(1 - alpha(i)^2 - beta(j)^2); 
        end
    end
end
% angular spectrum based formula
U1 = ifft2(ifftshift(fftshift(fft2(layer)).*exp(1i*k.*gamma_cust.*z))); 
I1 = (1/(16*pi)).*(U1.*conj(U1)); 
% show results 
figure; imshow(layer, []); title('circular aperture'); 
figure; imshow(I1, []); title('angular spectrum based diffraction');