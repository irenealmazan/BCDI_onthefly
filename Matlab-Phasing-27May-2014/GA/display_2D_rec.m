function  display_2D_rec(parent_dir,flip,flop)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    flip;
    if isempty(flip), flip=0;end
catch
    flip=0;
end

try
    flop;
    if isempty(flop), flop=0;end
catch
    flop=0
end

save_dir=parent_dir;

amp_f=rdir([parent_dir,'**/*AMP.rec']);
ph_f=rdir([parent_dir,'**/*PH.rec']);
pn=load_rec(amp_f.name,ph_f.name,flip);

if flop ~= 0,pn=flipdim(pn,2);end

nx=size(pn,2);
ny=size(pn,1);

lw=1.5;

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4)),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4)) );
imagesc(ph)
h=colorbar('location','EastOutside','fontsize',0,'fontweight','bold')
caxis([-pi, pi])
saveas(fh, [save_dir,'Ph-c'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Ph-c']);


fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
imagesc(abs(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4)))
saveas(fh, [save_dir,'Amp-c'],'epsc');
print(fh, '-dpng','-r300', [save_dir,'Amp-c']);




end

