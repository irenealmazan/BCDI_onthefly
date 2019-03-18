% this script display the result of the new support:

figure;

subplot(7,2,1);imagesc(abs(rho(:,:,30)));axis image
subplot(7,2,2);imagesc(angle(rho(:,:,30)));axis image
title('rho');

subplot(7,2,3);imagesc(abs(rho_fft(:,:,30)));axis image
subplot(7,2,4);imagesc(angle(rho_fft(:,:,30)));axis image
title('fft(rho)');

subplot(7,2,5);imagesc(abs(gauss_to_convolve(:,:,30)));axis image
subplot(7,2,6);imagesc(angle(gauss_to_convolve(:,:,30)));axis image
title('gaussian');

subplot(7,2,7);imagesc(abs(rho_fft_gauss(:,:,30)));axis image
subplot(7,2,8);imagesc(angle(rho_fft_gauss(:,:,30)));axis image
title('fft(rho)*gaussian');

fft_test = fftshift(ifftn(fftshift(rho_fft_gauss)));

subplot(7,2,9);imagesc(abs(fft_test(:,:,30)));axis image
subplot(7,2,10);imagesc(angle(fft_test(:,:,30)));axis image
title('ifft(fft(rho)*gaussian)');

result1 = rho.*fft_test;

subplot(7,2,11);imagesc(abs(result1(:,:,30)));axis image
subplot(7,2,12);imagesc(angle(result1(:,:,30)));axis image
title('rho*ifft(fft(rho)*gaussian)');

subplot(7,2,13);imagesc(abs(new_support(:,:,30)));axis image
subplot(7,2,14);imagesc(angle(new_support(:,:,30)));axis image