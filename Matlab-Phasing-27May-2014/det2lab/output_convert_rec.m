function  output_convert_rec( ddir ,conj_ref,zero_ph,sample_pixel,apply_sup,newxyz)
%jclark
%load the output amp and phase from the normal matlab phasing
%and do the cooirdinate transfrom, output some images and
%vtk files.  
%sign change is necessary for phase to get the directiosn correct
% i.e x-> y| and same for vtk.

%conj reflect the reconstruction? (=1 yes, =0 no)
if exist('conj_ref') ~= 1,conj_ref=0;end
if exist('zero_ph') ~=1,zero_ph=0;end
if exist('apply_sup') ~=1,apply_sup=0;end


% get the paramter file
pfile=rdir([ddir,'*PARAMS.mat']);


%get the name and dir
prefix=pfile.name(numel(ddir)+1:end-10);

%load the params and reconstruction
load(char(pfile.name))
pn=load_rec_from_dir(ddir);

if exist('sample_pixel') == 1,params.sample_pixel=sample_pixel;end

if exist('newxyz') ~=1,params.newxyz=[];else params.newxyz=newxyz;end

%load the support
try
    sup = load_sup_from_dir(ddir);
catch
    sup=ones([size(pn)]);
end
%cehck for conj ref
if conj_ref == 1,
    disp('Conjugating and reflecting Lab frame output....')
    pn=conj_reflect(pn);
    sup=conj_reflect(sup);
end

if zero_ph == 1
    disp('Zeroing phase....')
    pn=zero_phase(pn);
end

%

%
save_dir=ddir;

if isdir([save_dir,'/Lab/']) == 0,mkdir([save_dir,'/Lab/']);end

%do the transfrom and save vtk
disp('<<<<<<< Performing coordinate transform >>>>>>>')
disp(' ')
disp('Changing phase sign to be consistent with geometry....')
ci=complex(0,1);
pn=abs(pn).*exp(-ci*angle(pn));
disp(' ')
disp('Transforming object....')
[array params]=det2lab(pn,params);
array=array./abs(max(array(:)));
plot_amp_phase(array,[save_dir,'/Lab/'])    
disp('Done....')

disp('Transforming support....')
[sup params]=det2lab(sup,params);
plot_amp_phase(round(abs(sup)),[save_dir,'/Lab/S'])
disp('Done....')
disp('Saving output to .vtk file....')

if apply_sup == 1,array=array.*shrink_wrap(abs(array),.1,1);end

disp('Saving transformed array to .rec file....')
save([save_dir,prefix,'LAB.rec'],'array')  
savevtk2scalar(flip_all_dim(abs(array)),[save_dir,'Amp-Phase.vtk'],flip_all_dim(angle(array)),params.sample_pixel)    

array=sup;
save([save_dir,prefix,'LAB-S.rec'],'array')
disp('Done....')

savevtk2scalar(flip_all_dim(abs(sup)),[save_dir,'Support.vtk'],[],params.sample_pixel)    
disp('Done....')

if params.pcdi == 1
    disp('Transforming MCF....')
    pfile=rdir([ddir,'**-COH.rec']);
    load(char(pfile.name),'-mat');
    [coh params]=det2lab(array,params);
    disp('Done....')
    savevtk2scalar(flip_all_dim(abs(array)),[save_dir,'Amp-Phase-Coh.vtk'],flip_all_dim(angle(array)),params.sample_pixel)    
    array=coh;save([save_dir,prefix,'LAB-C.rec'],'array');
    plot_amp_phase(coh,[save_dir,'/Lab/Coh'])
end

%do the figs output
%if isdir([save_dir,'Lab/']) ==0,mkdir([save_dir,'Lab/']);end

%plot_amp_phase(array,[save_dir,'Lab/'])

try
    [ params_resn ] = determine_resolution(abs(array));
    params_resn.dsx=params.sample_pixel;params_resn.dsy=params.sample_pixel;params_resn.dsz=params.sample_pixel;
    output_resolution(params_resn,[save_dir,'Lab/'])
catch
    disp('ERROR calculating resolution....') 
end


clear support

close all

end

function plot_amp_phase(pn,save_dir)

nx=size(pn,2);
ny=size(pn,1);
nz=size(pn,3);

xx=[1,nx];%[nx/2-nx/3,nx/2+nx/3];
yy=[1,ny];%][ny/2-ny/3,ny/2+ny/3];
zz=[1,nz];%][nz/2-nz/3,nz/2+nz/3];

lw=1.5;
ph_range=[-pi,pi];
phase=atan2(imag(pn),real(pn) );
amp=abs(pn);

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(amp,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
%saveas(fh, [save_dir,'Amp-xy'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Amp-xy']);
%plot2svg([save_dir,'Amp-xy'],fh)
save_figure_mult(fh,save_dir,'Amp-xy')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'xz' );
imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
%saveas(fh, [save_dir,'Amp-xz'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Amp-xz']);
save_figure_mult(fh,save_dir,'Amp-xz')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(amp,'zy' );
imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
%caxis(ph_range)
%saveas(fh, [save_dir,'Amp-zy'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Amp-zy']);
save_figure_mult(fh,save_dir,'Amp-zy')

%****************

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
%ph=atan2(imag(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))),real(pn(ny/2-ny/4:ny/2+ny/4,nx/2-nx/4:nx/2+nx/4,round(nz/2))) );
ph = extract_3D_slice(phase,'xy' );
imagesc(ph(yy(1):yy(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
%saveas(fh, [save_dir,'Ph-xy'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Ph-xy']);
save_figure_mult(fh,save_dir,'Ph-xy')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(phase,'xz' );
imagesc(ph(zz(1):zz(2),xx(1):xx(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
%saveas(fh, [save_dir,'Ph-xz'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Ph-xz']);
save_figure_mult(fh,save_dir,'Ph-xz')

fh = figure ; % returns the handle to the figure object
set(fh, 'color', 'white'); % sets the color to white 
ph = extract_3D_slice(phase,'zy' );
imagesc(ph(yy(1):yy(2),zz(1):zz(2)));
h=colorbar('location','EastOutside','fontsize',20,'fontweight','bold');
caxis(ph_range);
%saveas(fh, [save_dir,'Ph-zy'],'epsc');
%print(fh, '-dpng','-r300', [save_dir,'Ph-zy']);
save_figure_mult(fh,save_dir,'Ph-zy')

end