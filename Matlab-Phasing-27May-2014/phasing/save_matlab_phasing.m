function  save_matlab_phasing(pnm,support,data,coh,chi,params,save_dir  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
disp('-----------------------------------------------------------')
disp('preparing reconstruction for saving....')
disp('')
disp('centering arrays....')
disp('')
%[pn xyz]=center_array(pnm);

axis_save=0;

try
    params.do_coord_tform;
catch
    params.do_coord_tform=1; %coord tform, =1 yes, =0 no
end

%%
if ndims(pnm) == 3
    
    %shift based on max value, required if it spans across the edge of
    %the array before COM
    [pnm yxz]=center_array(pnm);
    support=circshift(support,yxz);
    
    %do COM centering, will already be close
    xyz=center_of_mass(abs(pnm).*support);
    xyz=-1*round([xyz(2),xyz(1),xyz(3)]);
    pn=circshift(real(pnm),xyz)+1i*circshift(imag(pnm),xyz);
    support=circshift(support,xyz);
    %coh=circshift(coh,xyz);  %not sure wht this was here, don't need to do
    
    %set COM phase to zero, use as a reference
    disp('removing phase offset....')
    disp('')
    sz=size(pn);
    i=round(sz(1)/2);
    j=round(sz(2)/2);
    k=round(sz(3)/2);
    phi0=atan2(imag(pn(i,j,k)),real(pn(i,j,k)));
    pn=pn*exp(-1i*phi0);
    
    if params.recenter_rec == 1
        %do a first approximation to remove phase ramp
        disp('removing phase ramp....')
        disp('')
        pn = remove_ramp_pn_ups(pn,3);
%         amp=abs(pn);
%         damp=fftshift(fftn(fftshift(amp)));
%         [l,m,n]=ind2sub(size(damp),find( abs(damp) == max(max(max(abs(damp))))));
%         dpn=circshift(fftshift(fftn(fftshift(pn))),[l-i,m-j,n-k]);
%         pn=fftshift(ifftn(fftshift(dpn)));
    end
end
if ndims(pnm) == 2
    xyz=center_of_mass((abs(pnm).*support));
    xyz=-1*round([xyz(2),xyz(1)]);
    pn=circshift(real(pnm),xyz)+1i*circshift(imag(pnm),xyz);
    support=circshift(support,xyz);
    %coh=circshift(coh,xyz);

    %set COM phase to zero, use as a reference
    disp('removing phase offset....')
    disp('')
    sz=size(pn);
    i=round(sz(1)/2);
    j=round(sz(2)/2);
    phi0=atan2(imag(pn(i,j)),real(pn(i,j)));
    pn=pn*exp(-1i*phi0);

    %do a first approximation to remove phase ramp
%     disp('removing phase ramp....')
%     disp('')
%     amp=abs(pn);
%     damp=fftshift(fftn(fftshift(amp)));
%     [l,m]=ind2sub(size(damp),find( abs(damp) == (max(max(abs(damp))))));
%     dpn=circshift(fftshift(fftn(fftshift(pn))),[l-i,m-j]);
%     pn=fftshift(ifftn(fftshift(dpn)));
end
%%
if ndims(pnm) == 3
    %flip arrays for correct orientation on output
    try
        flip=params.flip;
    catch
        flip=1; %flip == 1 wil be correct
        params.flip=1;
    end
    
    try
       params.det_orient; 
    catch
       params.det_orient='xy';
    end

    if strcmp(params.det_orient,'yx')
        flip=0;end     %flips reconstruction up/down

    try 
        flop=params.flop;    
    catch
        flop=0;     %flips left/right
    end
    axis_save=1;

    if flip == 1
       nf=size(support);
       nf=nf(3);

       for ff = 1:nf, 
           support(:,:,ff)=flipud(support(:,:,ff));
           pn(:,:,ff)=flipud(pn(:,:,ff));
           data(:,:,ff)=flipud(data(:,:,ff));
           if sum(size(params.coh))>10
           coh(:,:,ff)=flipud(coh(:,:,ff));
           end
       end
    end
    if flop == 1
       nf=size(support);
       nf=nf(3);

  
       for ff = 1:nf, 
           support(:,:,ff)=flipud(rot90(support(:,:,ff),2));
           pn(:,:,ff)=flipud(rot90(pn(:,:,ff),2));
           data(:,:,ff)=flipud(rot90(data(:,:,ff),2));
           if sum(size(params.coh))>10
           coh(:,:,ff)=flipud(rot90(coh(:,:,ff),2));
           end
       end
           
    end
end
%%
files=params.files;
ALG1=params.ALG1;
ALG2=params.ALG2;
iterations=params.iterations;
pcdi=params.pcdi;
sw=params.sw;
seq_let=params.seq_let;
pcdi_type=params.pcdi_type;



%% save some things as matlab files.  despite the .rec, they are matlab.
%% allows them to be distinguished easily for transfroming to .vtk files
file = create_save_name( files,ALG1,ALG2,iterations,pcdi,sw,seq_let,pcdi_type,params.GPU,params );

disp(['save name - ',file])

save_dir=[save_dir,file];

if isdir(save_dir) ==0,mkdir(save_dir);end

save_dir=[save_dir,'/'];

%save amp
disp('saving reconstructed object amplitude....')
array=abs(pn)/max(max(max(abs(pn))));
%save([save_dir,file,'-AMP.mat'],'array')
disp('done...')

%save_phase
disp('saving reconstructed object phase....')
array=atan2(imag(pn),real(pn));
%save([save_dir,file,'-PH.mat'],'array')
disp('done...')
    
%save support
disp('saving support....')
array=support;
%save([save_dir,file,'-SUP.mat'],'array')
disp('done....')

%save data if save_data='YES'
data_save=0;
%if strcmp(params.save_data,'YES'), data_save=1;end
%if strcmp(params.save_data,'yes'), data_save=1;end

if data_save ==1,
    disp('saving data....')
    array=data;
 %   save([save_dir,file,'-DATA.mat'],'array')
    disp('done....')
else disp('not saving data....'),end

if axis_save ==1,
    disp('saving axis....')
    array=generate_axis(sz(2),sz(1),sz(3),20);
    %array=rescale_rec(array,params.det_px*params.binning(1),params.det_px*params.binning(2),params.dth,params.arm,params.lam);
  %  save([save_dir,file,'-AXIS.mat'],'array')
    disp('done....')
else disp('not saving axis....'),end
    
%save coh
if pcdi == 1,
    disp('saving recovered coherence function....')
    array=center_array(coh);
    save([save_dir,file,'-COH.mat'],'array'); 
    disp('done....')
else disp('pcdi off, not saving coherence function....'),end

%save error and params file
disp('saving chi...')
save([save_dir,file,'-ERROR.mat'],'chi')
disp('done....')
disp('saving parameter file....')
%save([save_dir,file,'-PARAMS.mat'],'params')
disp('done....')

%save a text file with the params as well
disp('outputting parameter .txt file....')
%reconstruction_output([save_dir,file,'-PARAMS.txt'],params)
disp('done....')
%%
disp('outputting pythonphasing.py for use with mayavi2....')
%phasing_params_out( params,save_dir );  %write a python file with the phasing 
disp('done....')                        %params for mayavi to read when
                                        %loading the reconstructions (.rec
                                        %files)
                                        
%%
disp('outputting MATLAB_to_vtk_Qv.py to generate .vtk files (inc Q vector) &')
disp('this will also output amp-phase.vtk containing both (amp & phase) in one file....')
%output_python_script(save_dir,'MATLAB_to_vtk_Qv.py');
%MATLAB_to_vtk_py(save_dir);            %outputs a python script that converts 
disp('done....')                        %matlab reconstruction to vtk
                                        %uses the ##PARAMS.mat file to load
                                        %experimental paramters.  these can
                                        %be overridden by manually changing
                                        %them in the outputted file,
                                        % MATLAB_to_vtk.py
%%                                        
disp('outputting Mayavi_PhaseIso.py to generate Mayavi images....')
%output_python_script(save_dir,'Mayavi_PhaseIso.py');
%MATLAB_to_vtk_py(save_dir);            %outputs a python script that converts 
disp('done....')                        %matlab reconstruction to vtk
                                        %uses the ##PARAMS.mat file to load
                                        %experimental paramters.  these can
                                        %be overridden by manually changing
                                        %them in the outputted file,
                                        % MATLAB_to_vtk.py
                                        %%                                        
disp('outputting Mayavi_AmpIso.py to generate Mayavi images....')
%output_python_script(save_dir,'Mayavi_AmpIso.py');
%MATLAB_to_vtk_py(save_dir);            %outputs a python script that converts 
disp('done....')                        %matlab reconstruction to vtk
                                        %uses the ##PARAMS.mat file to load
                                        %experimental paramters.  these can
                                        %be overridden by manually changing
                                        %them in the outputted file,
                                        % MATLAB_to_vtk.py
%%                                        
disp('outputting Mayavi_scene_script_Ver1-2.py to generate Mayavi images....')
%output_python_script(save_dir,'Mayavi_scene_script_Ver1-2.py');
%MATLAB_to_vtk_py(save_dir);            %outputs a python script that converts 
disp('done....')                        %matlab reconstruction to vtk
                                        %uses the ##PARAMS.mat file to load
                                        %experimental paramters.  these can
                                        %be overridden by manually changing
                                        %them in the outputted file,
                                        % MATLAB_to_vtk.py
%%
disp('outputting Make_moviePh.py to generate movie (requires ffmpeg)....')
%output_python_script(save_dir,'Make_moviePh.py');          
disp('done....') 
%%
disp('outputting Make_movieAmp.py to generate movie (requires ffmpeg)....')
%output_python_script(save_dir,'Make_movieAmp.py');          
disp('done....') 
%%
% %no need for combine.py, MATLAB_to_vtk_Qv.py does all this now.  
%disp('outputting Make_movie.py to generate movie (requires ffmpeg)....')
%output_python_script(save_dir,'combine.py');          
%disp('done....')                        

%% plot some stuff

if ndims(pnm) == 3,save_display_rec3D(pn,support,save_dir );end %save pictures
if ndims(pnm) == 2,save_display_rec2D(pn,support,save_dir );end

%do the coord transform
if ndims(pnm) == 3
    if sum([params.save_to_vtk,params.save_det2lab,params.calc_resolution]) > 0
        
        if isdir([save_dir,'/Lab']) ~= 1,mkdir([save_dir,'/Lab']);end
        
        if params.do_coord_tform == 1
            disp('<<<<<<< Performing coordinate transform >>>>>>>')
        else
            disp('<<<<<<< NOT Performing coordinate transform >>>>>>>')
        end    
        disp(' ')
        disp('Changing phase sign to be consistent with geometry....')
        ci=complex(0,1);
        pn=abs(pn).*exp(-ci*angle(pn));
                    
        if params.do_coord_tform == 1
            disp(' ')
            disp('Transforming object....')
            [array params]=det2lab(pn,params);
        else 
            array=pn;
        end
        
        array=array./(max(abs(array(:))));
        disp('Done....')
        plot_amp_phase(array,[save_dir,'/Lab/'])    
        
       
        if params.do_coord_tform == 1
            disp('Transforming support....')
            [support params]=det2lab(support,params);
            disp('Done....')
        end
        
        if params.pcdi == 1
            if params.do_coord_tform == 1
                disp('Transforming MCF....')
                [coh params]=det2lab(coh,params);
                disp('Done....')
            end
            plot_amp_phase(coh,[save_dir,'/Lab/Coh'])
        end
                        
        
        if params.calc_resolution == 1
            try
                [ params_resn ] = determine_resolution(abs(array));
                params_resn.dsx=params.sample_pixel;params_resn.dsy=params.sample_pixel;params_resn.dsz=params.sample_pixel;
                output_resolution(params_resn,[save_dir,'/Lab/'])
            catch
                disp('ERROR calculating resolution....')
            end
        end
        if params.save_to_vtk == 1
            array=flip_all_dim(array);
            support=flip_all_dim(support);
            coh=flip_all_dim(coh);
            disp('Saving output to .vtk file....')
           % savevtk2scalar(abs(array),[save_dir,'Amp-Phase.vtk'],angle(array),params.sample_pixel)    
           % savevtk2scalar(abs(support),[save_dir,'Support.vtk'],[],params.sample_pixel)    
          %  if params.pcdi == 1,savevtk2scalar(abs(coh),[save_dir,'Amp-Phase-Coh.vtk'],angle(coh),params.sample_pixel);end 
            disp('Done....')
        end
        if params.save_det2lab == 1
            disp('Saving transformed array to .rec file....')
            save([save_dir,file,'-LAB.mat'],'array')  
            array=support;
          %  save([save_dir,file,'-LAB-S.mat'],'array')
          %  if params.pcdi == 1,array=coh;save([save_dir,file,'-LAB-C.rec'],'array');end
            disp('Done....')
            clear support
            clear coh
        end
        
        disp('<<<<<<< Finished coordinate transform >>>>>>>')
    end
end


disp('-----------------------------------------------------------')
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

%****************

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

end
