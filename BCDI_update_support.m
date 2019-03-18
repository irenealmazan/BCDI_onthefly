function [new_support] = BCDI_update_support(rho)
    % this function updates the support by convolving the updated rho
    % with a gaussian function and then selecting the points which have
    % an intensity above a threshold
    global X Y Z
    
    [Npix_x,Npix_y,Npix_z] = size(X);
    
    threshold = 0.20; % 20% of the signal
    
    rho_fft = fftshift(fftn(fftshift(rho)));

    FWHM_x = 0.1;
    FWHM_y = 0.1;
    FWHM_z = 0.1;
    
    gauss_to_convolve = exp(-FWHM_x.*X.^2/Npix_x-FWHM_y.*Y.^2/Npix_y-FWHM_z.*Z.^2/Npix_z);
    
    rho_fft_gauss = rho_fft.*gauss_to_convolve; 
    
    rho_conv_gauss = fftshift(ifftn(fftshift(rho_fft_gauss)));
    
    max_of_int = max(max(max(rho_conv_gauss)));
    
    nonzeros_elements = find(rho_conv_gauss> threshold*max_of_int);

    new_support = zeros(Npix_x,Npix_y,Npix_z);
    new_support(nonzeros_elements) = 1;

end