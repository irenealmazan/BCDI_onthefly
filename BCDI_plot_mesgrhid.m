% plot the meshgrid which has been created:

[X_square,Y_square,Z_square] = meshgrid([-Npix_x/2 Npix_x/2-1]*d2_bragg,[-Npix_y/2 Npix_y/2-1]*d2_bragg,[-Npix_z/2 Npix_z/2-1]*d2_bragg);
X_square_toplot = X_square(:);
Y_square_toplot = Y_square(:);
Z_square_toplot = Z_square(:);

figure(3);
hold on;
scatter3(X_square_toplot,Y_square_toplot,Z_square_toplot)
%}