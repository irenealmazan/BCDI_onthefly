function sz = plot_amp_phase_v3(pn,save_dir,ph_range,sz,params)

if ischar(save_dir) == 0
   
    fdir=strtrim(char(save_dir(1)));
    sdir=strtrim(char(save_dir(2)));
    tdir=strtrim(char(save_dir(3)));
    
else
    
    fdir=save_dir;
    sdir=save_dir;
    tdir=save_dir;
    
end

[fdir0 fout]=extract_dir_from_string(fdir);
[sdir0 sout]=extract_dir_from_string(sdir);
[tdir0 tout]=extract_dir_from_string(tdir);

%check directorys exist
make_directory({sdir0,fdir0,tdir0});
%

nx=size(pn,2);
ny=size(pn,1);
nz=size(pn,3);

switch sum(sz(:))
    
    case 0
        fact=1.5;
        num_x=abs(pn(round(ny/2),:,round(nz/2)));
        num_y=abs(pn(:,round(nx/2),round(nz/2)));
        num_z=abs(pn(round(ny/2),round(nx/2),:));
        ind_x=(num_x >= .1*max(num_x));
        ind_y=(num_y >= .1*max(num_y));
        ind_z=(num_z >= .1*max(num_z));
        sx=fact*sum(ind_x(:));
        sy=fact*sum(ind_y(:));
        sz=fact*sum(ind_z(:));

        xx=round([nx/2-sx/2,nx/2+sx/2]);
        xx(xx < 1)=1;
        xx(xx > nx)=nx;
        yy=round([ny/2-sy/2,ny/2+sy/2]);
        yy(yy < 1)=1;
        yy(yy > ny)=ny;
        zz=round([nz/2-sz/2,nz/2+sz/2]);
        zz(zz < 1)=1;
        zz(zz > nz)=nz;
        sz=[xx,yy,zz];
    case 1
        xx=[1,nx];%[nx/2-nx/3,nx/2+nx/3];
        yy=[1,ny];%[ny/2-ny/3,ny/2+ny/3];
        zz=[1,nz];%[nz/2-nz/3,nz/2+nz/3];
        sz=[xx,yy,zz];
    otherwise
        xx=[sz(1),sz(2)];
        yy=[sz(3),sz(4)];
        zz=[sz(5),sz(6)];
end
    


lw=1.5;

phase=atan2(imag(pn),real(pn) );
amp=abs(pn);



fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(amp,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
axis equal
axis tight
%caxis(ph_range)

box off

try
    params.ampmap;
    colormap(bone);
    %colormap(flipud(colormap));
    box off
    axis off
end

set(gca, 'LooseInset', get(gca, 'TightInset'));

%saveas(fh, [fdir,'Amp-xy'],'epsc');
print(fh, '-dpng','-r300', [fdir,'Amp-xy']);


    
if isreal(pn) ~= 1
    colormap(jet)
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    %ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
    ph = extract_3D_slice(phase,'xy' );
    imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
    h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    axis tight
    
    if isempty(ph_range) ~= 1
        caxis(ph_range);
    end
    
    box off
    %axis off
    
    %saveas(fh, [fdir,'Ph-xy'],'epsc');
    print(fh, '-dpng','-r300', [fdir,'Ph-xy']);
    
    colormap(jet)
    fh = figure ; % returns the handle to the figure object
    set(fh, 'color', 'white'); % sets the color to white 
    %ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
    pn = (extract_3D_slice(pn,'xy' ));
    imagesc(c2image(pn(yy(1):yy(2),xx(1):xx(2))));
    %h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
    axis equal
    axis tight
    %caxis(ph_range);
    
    box off
    %axis off
    
    %saveas(fh, [fdir,'Ph-xy'],'epsc');
    print(fh, '-dpng','-r300', [fdir,'xy']);

    
end
% fh = figure ; % returns the handle to the figure object
% set(fh, 'color', 'white'); % sets the color to white 
% imagesc(abs(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))))
% h=colorbar('location','EastOutside','fontsize',0,'fontweight','bold');
% saveas(fh, [save_dir,'-Amp-c'],'epsc');
% print(fh, '-dpng','-r300', [save_dir,'-Amp-c']);


end

