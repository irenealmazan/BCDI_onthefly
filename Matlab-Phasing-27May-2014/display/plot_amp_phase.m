function plot_amp_phase(pn,save_dir,ph_range)

nx=size(pn,2);
ny=size(pn,1);
nz=size(pn,3);

xx=[nx/2-nx/3,nx/2+nx/3];
yy=[ny/2-ny/3,ny/2+ny/3];
zz=[nz/2-nz/3,nz/2+nz/3];


lw=1.5;

phase=atan2(imag(pn),real(pn) );
amp=abs(pn);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(amp,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
saveas(fh, [save_dir,'Amp-xy'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Amp-xy']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'xz' );
imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
saveas(fh, [save_dir,'Amp-xz'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Amp-xz']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'zy' );
imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
saveas(fh, [save_dir,'Amp-zy'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Amp-zy']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(phase,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
saveas(fh, [save_dir,'Ph-xy'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Ph-xy']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(phase,'xz' );
imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
saveas(fh, [save_dir,'Ph-xz'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Ph-xz']);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(phase,'zy' );
imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
saveas(fh, [save_dir,'Ph-zy'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Ph-zy']);

% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% imagesc(abs(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))))
% h=colorbar('location','EastOutside','fontsize',0,'fontweight','bold');
% saveas(fh, [save_dir,'-Amp-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Amp-c']);


end

