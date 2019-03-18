%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlabphasing_ver1.1.m 
%2D/3D phase retrieval - Jesse Clark, LCN, UCL October-November 2010
%                       jesse.clark@ucl.ac.uk, jessenclark@gmail.com
clc
clear all
params.version='Matlab phasing version 1.1 - Nov 2013';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

snums= [46];


for tt = 1

for tq = 2:5


%% Determining the current directory
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data preperation
params.data_dir=dir;            %assumes data is in the same dir as this file
                                %otherwise specify the data directory.
                                % at the moment it will create a save dir with a name that is 
                                % generated (see bottom of page) according to phasing params.  it will be
                                % saved in the directory that this phasing file is in.  this can be over-
                                % ridden by simply specifiying another save_dir.  
params.save_dir=dir;

                                % the data files and background files to load.  if no bg files are needed
                                % put back={}.  Accepts .spe files and also matlab save files .mat.  if
                                % matlab is used it is assumed the array is 3d.  all other options work the
                                % same.

%params.files={'oleg1113_S0566.tif','oleg1113_S0567.tif','oleg1113_S0568.tif'};%
%params.files={'S0400.mat','S0401.mat','S0402.mat'};
params.files={['S0' int2str(snums(tt)) '.tif']};
params.back={};%if no bg required, leave empty {}.  don't use {''}.

params.itnum = 1;               %change this for each recon of same data set
 
params.seq_let=['r' num2str(tq)]; %sequence letter, second part of save name

%params.comments2 = '40H'; %first part of save name

params.comments2 = ['gold' num2str(snums(tt))];

params.comments='';      %write comments here and this will be saved in the params file

params.binning=[1,1];           %binning for x and y

params.skipping=0 ;             %will bin (=0), will skip points (=1)

params.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=250;            %min data threshold.  below this is removed.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.
                                
params.schot_th=0;           %secondary threshold applied AFTER binning
params.subtract_dc=0;           %leave this as 0 (unless you know what you are doing)
params.no_center=0;             %no centering of data -> no_center=1
params.no_fft_pad=1;            %no padding for fft -> no_fft_pad=1
params.pad_ptych=0;             %leave this as 0 (unless you know what you are doing)
params.no_hist=1;               %plot the histograms of the data (=0)
params.bg_mult=1;               %multiplication to apply to the bg file (if exp time is different)
params.save_data='YES';


params.data_only='NO';        
                                %will save or exit after data prep. ='NO' no save and no exit,         
                                %'YES' will save and continue, 'YES-EXIT'
                                %will save and then stop the script

params.nnc=[0,0,0,0,0,0];
                                % initial cropping of data before binning.
                                %eg. nnc=[0,0,-10,-10,5,5] will do nothing to x,
                                %will crop 10 pixels off each end in y and will pad
                                %5 pixels to each end in z.  set nnc=0 to do
                                %nothing.
params.do_2D=0;                 %will take the central slice and do 2D phase retrieval if =1                                
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% No inputs required here.  Doing the data preperation.                     
[data params] = load_MP_data(params);   %load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Experimental Geometry
% used for displying the reconstructions

params.det_px=22.5e-6;               %detector pixel size (m)
params.lam=.1358;                   %wavelength (nm)
params.delta=32.6;                  %delta (degrees)
params.gam=10.4;                     %gamma (degrees)
params.arm=1.5;                      %camera distance (m)
params.dth=(.15+.15)/30;                      %angular step size
params.dtilt=0;


params.spx=params.arm*params.lam/params.nn(1)/params.binning(1)/params.det_px;   %sample pixel size
params.spy=params.arm*params.lam/params.nn(2)/params.binning(2)/params.det_px;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Algorithm selection and control      
params.iterations= 821;             % Iterations to perform.

params.start_guess='random-data';          
                                    % start guess for iterative procedure. 'flat' uses the support
                                    % 'random' uses a random object the
                                    % size of the support.  'random-data'
                                    % uses the data with a random phase.
                                    % if the full path to a previous
                                    % reconstruction is supplied it will
                                    % use that.  it will align the
                                    % reconstruction with the data as well.
                                      
                                    %select one or two algorithms.  Choices are 'ER','SF','HIO','DM',
                                    %'RAAR','HPR','ASR'.  Switching between algorithms is determined
                                    %by trigger. to use one algorithm set, ALG2=ALG1
params.ALG1='ER' ;                  %Algorithm 1.  must be in capitals
params.ALG2='HIO' ;                 %Algorithm 2.  must be in capitals

params.trigger=[10,40,50,80,90,120,130,160,170,200,210,240,250,280,290,320,330,360,370,400,410,440,450,480,490,520,...
	530,560,570,600,610,640,650,670,680,690,710,720,730,750,760,770,780,800];

%params.trigger=5:5:800;
                                    %trigger is an array of iteration numbers 
                                    %that determine when to switch
                                    %algorithms.  ie. trigger=[10,100,110,150] will do
                                    %ALG1 until iteration 10 then switch to
                                    %ALG2 until iteration 100 then it will
                                    %switch back to ALG1 at 110 etc etc. until all iterations
                                    %are complete.  There is no limit to how many times it
                                    %can change.  The number of iterations is set by
                                    %iterations though
                
params.beta=.9;                     %HIO,RAAR,HPR alg parameter
params.twin=1;                      %set this to do another support (twin=1) that cuts out half the array
                                    %it gets rid of the twin.  will only do it once at an
                                    %iteration specified in the algorithm by twin_it.
params.twin_it=15;                   %set twin=0 to turn off       


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Support parameters, including shrinkwrap and phase constraint
% Parameters controlling the support.  To turn shrinkwrap OFF either set
% threshold=0 OR set the elements of shrink_trigger <0.


%0.156 is 40 pixels
%max is 40 pixels in x
%max is 25 pixels in y
%sx=40;
sx=0.25*params.nn(1);
sy=0.25*params.nn(2);   
sz=0.3*params.nn(3);   %x,y,z support size (pixels)
%sy=40;%0.156*params.nn(2);                  %i like to set them relative to the array dimensions   
%sz=20;%0.30*params.nn(3);                  %but you can put in pixel values as well
 
params.threshold=0.1;                 %Shrinkwrap threshold (0<threshold<1).  threshold=0 will turn SW OFF.
params.sigma=1; %jesse says close to 1                        %Sigma of gauss kernal for smoothing in SW
params.sw_type='gauss';
                                       %'gauss' - uses a gaussian Kernal to smooth
                                       %'box' - uses a box kernal for smothing. 
                                       %'percent' - adjusts the threshold so the support is
                                       %            a constant number of voxels.  ie. .1 will make sure
                                       %            the support is 10% of the total voxels
                                       %'gauss_fill' - is the same as gauss but fills in any
                                       %               holes, basically makes the support convex
                                       %'gauss_percent' - same as percent
                                       %                  but smooths first
                                       %'percent-auto'- estimates the size
                                       %              from the
                                       %              autocorrelation and
                                       %              keeps this fixed for
                                       %              the support
                                       %'gauss-minarea' - same as gauss but
                                       %                will switch to 'gauss-percent' if
                                      %                the support volume falls below params.sw_min_area

params.sw_min_area=.02;                %for use with 'gauss-minarea'
                                       %as a fraction of the total, i.e .02
                                       %is 2 percent
                                       
                              
params.shrink_trigger=generate_trigger(5,4,805);
                                       %the same as trigger but for controlling 
                                       %shrinkwrap.  make all <0 to keep shrinkwrap 
                                       %off or set threshold=0.
                                
params.phase_range=[-pi/2,pi/2];       %Pixels with a phase outside this range will have there
                                       %amplitude->0. Done within the
                                       %initial support region or
                                       %shrinkwrap support region
                                       
params.phase_trigger=[1,40];         %same rules as for trigger and shrink_trigger
                                       %i.e phase_trigger=[1,70] truns it on at
                                       %iteration 1 and off at 70.
                                       %phase_trigger=[-1,-70] won't turn on at
                                       %all.
                                       
                                       
%B paramsters controlling if a low to high resolutionreconstruction will
%take place or if the shrinkwrap paramters change during iteration
%to enable put a 'lr' on the end of each algorithm choice, e.g. ERlr or
%HIOlr will enable this for tose algorithms.  these fields will be set to a
%defualt set of values if left out but the 'lr' must be present.
                                            %change the sw sigma value?
                                            %must be a array of values with
                                            %the an element for each
                                            %iteration
%params.nonGA_lres_sig=make_vals(params.iterations-100,params.sigma,3);
                                            %the st dev of a guasian mask 
                                            %for the data to reconstruct a
                                            %low to high res version, very
                                            %useful.  leave as >=1 to do
                                            %nothing
%params.nonGA_lres_det=make_vals(params.iterations-100,1,.1);
                                       
                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Partial coherence correction parameters                                 
%A, what things to do, x,y,z etc.  angle is only used for the 'gauss_#'
%corrections.  

params.pcdi=0;                   %do PCDI correction, 1 is ON, 0 is OFF
                                 
params.x_not_do=0;               %set = 1 to not do x
params.y_not_do=0;               %set = 1 to not do y
params.z_not_do=0;               %set = 1 to not do z
                                 
params.roi=[16,16,16];           %ROI pixel values to do pcdi min over

%B, algorithm type and options

params.pcdi_type='gauss_ps';         %type of pcdi correction. 
                                 %'gauss','gauss_ps' - pattern search, 'gauss_cm' constrained match - assume gauss coh func  
                                 %'lucy' - arbitarily shaped coh func.

params.use_k_1=0;                %uses I(k-1) to determine kernal, used to ensure uniqueness
params.use_2k_1=1;               %uses 2*I(k)-I(k-1) to determine kernal, used to ensure uniqueness, USE THIS
params.conv_kern_size=-1;        %kernal size for convolution, must be smaller than roi.  set =-1 to use same size as roi.                                
params.lucy_its=30;              %RL sub-iterations
params.use_fftconv=1;            %use fft based convolution (=1)
params.use_previous_coh=1;       %start each RL sub-iteration from the previous output (=1)
params.symmetrize_kernal=0;      %enforce symetry of kernal in det plane (will be enforced in sample since it is real)
                                 
params.update_its=15;               %every X iterations will update the kernal, saves time
params.start_pcdi=80;               %start doing PCDI correction after this many iterations
params.stop_pcdi=800; %when to stop the updating of the kernal
params.pcdi_norm=1;                 %normalise before doing minimisation, 1=yes,0=no
params.pcdi_ER=0;                   %do minimisation on FT[ PiS[pn] ]

%C, paramters if using 'gauss_ps','gauss_cm' or 'gauss'.  max
%sigmas,rotation etc.  don't worry if using 'lucy'

params.kernalxy=[0.5,0.5,0.5];           %guess values for intial kernal, sigma [x,y,z].  not that important
                                        %just don't make it too large.  ~.5 for
                                        %each is a good start guess

params.angle_not_do=0;                  %set = 1 to not do angle (used in 2D and 'gauss'
                                        %correction only.  rotates in the x-z
                                        %plane
params.angle=0;                         %initial guess at angle (x-z plane)
params.angle_start=0;                   %iteration to start

params.max_sigmas=[1.5,1.5,1.5];        %upper/lower limit for the gauss kernal sigma size
params.min_sigmas=[0,0,0.1];            %used in the minimisation.  only used when
                                        %pcdi_type is any schemes with 'gauss..'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. Guided paramaters
params.GA=0;                            %turns it on (=1)
params.GA_switch=0;                     %switch to breed_mode2 if not updated
params.generations=5;                   %do this many generations  
params.population=5;                    %with this many individuals

params.GA_metric='chi';           %how to decide what is the 'best'
                                        %us used for the breeding only
                                        %'chi' - usual chi squared (min
                                        %chi)
                                        %'sharpness' - {|p|^4 (min
                                        %sharpest)
                                        %'area' -  the area of the support (largest) 
                                        %'TV' - min of the TV norm
params.GA_metric_return='chi';    % how to decide what to return at the end
                                            
                                        
params.GA_threshold=0.0;                %to do a support after each generation rather than during
                                        %0=off, 0<th<1. 
params.GA_sig=.5;                       %size of the blur for the support after each generation
                                        
params.breed_mode1='sqrt_ab';           %primary breeding mode, 'none' will do no breeding
                                        
                                        
params.breed_mode2='none';              %mode to switch to if not updated

params.GA_lres=1;                       %start with low res data (=1) or always use full res (=0)
params.GA_lres_init=.1;                 %initial sig, 0-1 (array size=1)
params.GA_lres_genstop=4;               %at which generation should the full res start being used
params.GA_lres_pow=1;                   %how the sig scales with generation (1 - linear, 2 - quad etc)
params.GA_sig_max=3;                    %max support sig value, will decrease linearly until params.sigma is reached

params.GA_sig_custom=make_vals(params.generations,1.0,params.GA_sig_max);
                                        %custom sig values, leave empty,[], for automatic mode
params.GA_lres_custom=make_vals(params.generations,.5,params.GA_lres_init);     

params.GA_iterations= params.iterations;              %Iterations to perform.  this is used instead of the previously
                                        %specified iterations
                                        
params.GA_trigger= params.trigger;              %trigger for switching algorithms.  same as specified earlier
                                        %but if GA=1 will use this trigger.
                                        
params.GA_return='best';                %what to return, 'best' or 'avg' or 'avg-half'
                                        %'best' returns the best accroding
                                        %to the metric selected
                                        % 'avg' and 'avg-half' returns the
                                        % average using all or half the
                                        % population
                                        % can also use 'best-cull' or
                                        % 'avg-half-cull' which will remove
                                        % members of the population after
                                        % each generation see npop_stop and
                                        % GA_cull_remove
                                        % also can do 'avg-#' where # is
                                        % any number.
                                        
params.npop_stop=2;                     % minimum population size after culling
params.GA_cull_remove=2;                % how many to remove after each generation
                                        
params.generations=3;                   %do this many generations  
params.population=15;                                      
params.GA_random_seed=1;                %seed for random number generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. Miscellaneous pramaters
params.mod_const_pow=1;          %the power to use for the modulus constraint
                                 %ie sqrt(Im/Ik)^l, l=1 is default.
params.regularized_amp='none';   %use a statistically regularized amplitude constraint?
                                 %'none', 'gauss' , 'poisson','unifrom'.
                                 %set to 'none' if unsure.
params.silent=0;                 %suppresses output to the screen         
params.no_zero=0;                %let zero measured intensity values float, 0=no,1=yes
                                 
params.GPU=0;                    %triggers GPU acceleration, won't do anything
                                 %if jacket is not installed.  not really
                                 %supported anymore

params.save_every_iterate =0;    %will save every iterate.  to be used for
                                 %for making a movie of the reconstruction

params.iterate_avg_save=0;       %save the average iterate
params.iterate_avg_start=250;    %when to start averaging iterates
params.iterate_avg_its=50;        %number of iterates to average
params.iterate_avg_inter=1;      %interval between adding averages

params.save_to_vtk=0;            %save vtk in lab coords (yes=1)
params.save_det2lab=1;           %save matlab in lab coords (yes=1). Makes x,y transverse to beam direction
params.calc_resolution=0;        %calc resolution of final reconstruction (yes=1)

% END OF USER INPUTS % END OF USER INPUTS % END OF USER INPUTS % END OF USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other stuff.  NO need to change anything.  makes a parameter structure
%% to save all the necessary info.  is saved at the end of the
%% reconstruction.
% used for storing all the parameters to output into a text file and some
% intial algorithm parameters

params.sz=size(data);
params.size_init=[sx,sy,sz];
params.orig=params.pcdi;
params.do=[params.x_not_do,params.y_not_do,params.z_not_do];
params.coh=gauss_3D(11,11,11,params.kernalxy(1),params.kernalxy(2),params.kernalxy(3));
params.update=1;                 %tells the alg weather to update the coh 
                                 %func, no need to change this
ss=round([sy,sx,sz]);            %its row column major so y is before x    

[pn support] = make_intial_support(data,ss,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The phasing algorithm part.  itertive phasing handles the most part.
%% it calls modulus_projector and next_iterate to to the phasing.
%% modulus_projector contains all the pcdi stuff via determine_coh
[pnm support params_out chi DM_error coh] = iterative_phasing_general(pn,data,support, params);

%% Save the matlab reconstruction and output matlab files to be read
%% by mayavi.  outputs a phasingparams.py file and Matlab_to_vtk.py file

save_matlab_phasing(pnm,support,data,coh,chi,params_out,params.save_dir);

%% Display some things

display_rec3D(pnm,support,chi,DM_error);
%clear
close all;

%% END OF FILE %% END OF FILE 

end

beep

end

%%
%particle loading procedure
clear all; clc;

for snum = 879;


numscans=5;
%some scans are 64 and some are 48 so be aware
%iterates = zeros(256,256,48,numscans);
%cohers = zeros(size(iterates));
chis = zeros(1,numscans);
linup = 1; %try to shift to match


%do the loading into a 4D matrix
for i=1:numscans

datdir = ['/Users/a4894z/Desktop/recon/part4' int2str(snum) 'r' int2str(i) '/'];

fname = ['part4' int2str(snum) 'r' int2str(i) '-LAB.mat']; %phase and amp

fname2 = ['part4' int2str(snum) 'r' int2str(i) '-LAB-S.mat']; %support

%fname3 = ['part4' int2str(snum) 'r' int2str(i) '-LAB-C.mat']; %coherence

fname4 = ['part4' int2str(snum) 'r' int2str(i) '-ERROR.mat']; %error, just take final value


anew = load([datdir fname]); %load amp/phase
%cnew= load([datdir fname3]); %load coh
snew = load([datdir fname2]);
enew = load([datdir fname4]);


%check and see how consistent they are
[xlim ylim zlim] = size(anew.array);
figure;
imagesc(abs(anew.array(:,:,zlim/2))); zoom(4);

%isosurface(abs(anew.array),0.3)

%store them in the arrays

iterates(:,:,:,i) = anew.array;

%cohers(:,:,:,i) = cnew.array;

supports(:,:,:,i) = snew.array;

chitemp = enew.chi;

chis(i) = chitemp(end);


end

[ylim xlim zlim] = size(anew.array);


%now do the averaging
ref = iterates(:,:,:,1);

for qq=2:numscans

   
    cnj_rf=is_conj_ref(ref,iterates(:,:,:,qq))
            
    if cnj_rf ==1
        iterates(:,:,:,qq)=conj_reflect(iterates(:,:,:,qq));
    end 
    
    %try to line them up using essentially cross correlation
    if linup == 1
    [h k l] = register_3d_reconstruction(abs(ref),abs(iterates(:,:,:,qq)));
    next_c = sub_pixel_shift(iterates(:,:,:,qq),h,k,l);
    iterates(:,:,:,qq) = next_c;
    end
    
end

close all;

%exclude some values
%iterates(:,:,:,2)=[];

asum = mean(iterates,4); %average amp and phase at same time
%csum = mean(cohers,4);
ssum = mean(supports,4);
chit = mean(chis);
%put all these bad boys into a structure

TOT.asum = asum;
%TOT.csum = csum;
TOT.ssum = ssum;
TOT.chit = chit;

%save the average in a mat file
v43 = genvarname(['part4' int2str(snum)]);
eval([v43 '= TOT;']);
save(['part4' int2str(snum)],v43); %saves it as 40A.mat, x40A is asum

clear TOT asum csum ssum chit iterates cohers supports;

end

%%
%assuming you've loaded them all in and they're all in structures. Need to
%do this to get the iterates matrix which I use later

i=1;

for snums = [879];

v43 = genvarname(['part4' int2str(snums)]);
eval(['TOT=' v43;]);
[xlim ylim zlim] = size(TOT.asum);


%store all scans in iterates again, need for aligning
iterates(:,:,:,i) = TOT.asum;


% subplot(6,6,i); 
% %figure;
% imagesc(angle(TOT.asum(:,:,zlim/2)).*(abs(TOT.asum(:,:,zlim/2))>=0.25));
% zoom(4);


phsf =  2*pi*sqrt(3)/8.1; %in Angstroms now
pval2=0.03; %change between .03 and 0.2 depending on isosurf vs imagesc
phasesp = (angle(TOT.asum).*(abs(TOT.asum)>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,350); %use spacing of 20nm or 200 angstrom
clear uy uz;
uxnans = ux;

%normal slice method

% subplot(6,4,i); 
% %figure;
% imagesc(uxnans);
% caxis([-2.5e-4 2.5e-4]);
% zoom(4.5);

%contour slice method

% h = contourslice(abs(TOT.asum), [], [], [zlim/2-2 zlim/2-1 zlim/2 zlim/2+1 zlim/2+2],8);
% camva(24); camproj perspective;
% campos([-3, -15, 5]);
% set(gcf, 'Color', [.3, .3, .3], 'renderer', 'zbuffer');
% set(gca, 'Color', 'black', 'XColor', 'black', ...
% 'YColor', 'black', 'ZColor', 'black');
% box off;
% view(-38,16)
% 
% %caxis([-4e-4 4e-4]);
% zoom(3.5);
% axis off;
% %set(gca,'color','k')
% set(h, 'LineWidth', 2);
% set(gcf,'color','w')

%projection onto isosurface to see evolution on surface/core

% ux(isnan(ux)==1) = 0;
% [x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
% 
% figure;
% cdata = smooth3(ux,'box',[3 3 1]);
% p = patch(isosurface(x,y,z,asum,0.1));
% isonormals(x,y,z,asum,p);
% isocolors(x,y,z,cdata,p);
% set(p,'FaceColor','interp','EdgeColor','none')
% %daspect([1 1 1]);
% axis tight
% camlight(230,250); 
% lighting GOURAUD; 
% %title('4.7 1st Charge ','FontSize',16)
% colorbar;
% %set(gca, 'ydir','reverse')
% view(3); 
% axis off;


%figure;
%isosurface(abs(TOT.asum),0.25)

figure; 
%imagesc(angle(TOT.asum(:,:,zlim/2)).*(avgamp(:,:,zlim/2)>=0.25));

imagesc(angle(TOT.asum(:,:,zlim/2)).*(abs(TOT.asum(:,:,zlim/2))>=0.2));
zoom(4);

i=i+1;

end

asum = mean(abs(iterates),4);

%%
%if you want to line them up

ref = iterates(:,:,:,3);
linup=1;


for qq=1:2      
        
   
    cnj_rf=is_conj_ref(ref,iterates(:,:,:,qq))
            
    if cnj_rf ==1
        iterates(:,:,:,qq)=conj_reflect(iterates(:,:,:,qq));
    end 
    
    %try to line them up using essentially cross correlation
    if linup == 1
    [h k l] = register_3d_reconstruction(abs(ref),abs(iterates(:,:,:,qq)));
    next_c = sub_pixel_shift(iterates(:,:,:,qq),h,k,l);
    iterates(:,:,:,qq) = next_c;
    end
    
end




%%
%look at the properties of the average

asum = mean(abs(iterates),4); %see what the average amplitude does

for z=1:zlim
	%imagesc(asum(:,:,z))
	imagesc(angle(iterates(:,:,z,1)).*(asum(:,:,z)>=0.25))
	zoom(5)
	pause
	
end

%check the strain
for q=1:8
phsf =  2*pi*sqrt(3)/8.1; %in Angstroms now
pval2=0.25;
phasesp = (angle(iterates(:,:,:,q)).*(abs(iterates(:,:,:,q))>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,350); %use spacing of 20nm or 200 angstrom
clear uy uz;
uxnans = ux;

figure; imagesc(ux(:,:,zlim/2)); zoom(5);

end


for z=1:zlim
	imagesc(ux(:,:,z))
	zoom(5)
	pause
end


%check the phase of a particular reconstruction cut by the average amplitude

for z=1:zlim
	%imagesc(avgamp(:,:,z))
	imagesc(angle(iterates(:,:,z,3)).*(asum(:,:,z)>=0.2))
	pause

	
end

%%
%compute some quantitative things

fDC = [400 401 403 404 421:425 429 430 434 439 440 442 445 447:453];


for q=1:length(fDC)

clear phasesp ampp ux;
%strain map, was using 0.1 before
pval=0.2;
aeff2 = 8.14; %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

ampp = abs(iterates(:,:,:,q));
[xmax ymax zmax] = size(ampp);
phasesp = (angle(iterates(:,:,:,q)).*(ampp>=pval))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan; %only want gradients where you have data
[ux uy uz] = gradient(phasesp,160,160,370); %use spacing of 20nm or 200 angstrom
clear uy uz;


%way that gives interpretable quant metrics
%center is between 0.53 and 1
uxcent = (ux.*(ampp>=0.53));
uxcent(uxcent==0)=nan;

%check number of elements in your calculation
uxcent2 = uxcent;
uxcent2(isnan(uxcent2)==1) = 0; nnz(uxcent2);
clear uxcent2;


%edge is between .2 and 0.25
uxedge = (ux.*(ampp<=0.25));
uxedge(uxedge==0)=nan;

uxedge2 = uxedge;
uxedge2(isnan(uxedge2)==1) = 0; nnz(uxedge2);
clear uxedge2;

%calculate the quant metrics
uxqtot(:,q) = reshape(ux,xmax*ymax*zmax,1);
uxqcent(:,q) = reshape(uxcent,xmax*ymax*zmax,1);
uxqedge(:,q) = reshape(uxedge,xmax*ymax*zmax,1);


uxnmean(q) = nanmean(uxqtot(:,q));
uxnmeancent(q) = nanmean(uxqcent(:,q));
uxnmeanedge(q) = nanmean(uxqedge(:,q));

%average of abs
uxavgtot(q) = nanmean(abs(uxqtot(:,q)));
uxavgcent(q) = nanmean(abs(uxqcent(:,q)));
uxavgedge(q) = nanmean(abs(uxqedge(:,q)));
%standard deviation
uxstdtot(q) = nanstd(uxqtot(:,q));
uxstdcent(q) = nanstd(uxqcent(:,q));
uxstdedge(q) = nanstd(uxqedge(:,q));
%maximums
uxmaxtot(q) = max(abs(uxqtot(:,q))); %this also might be useful
uxmaxcent(q) = max(abs(uxqcent(:,q)));
uxmaxedge(q) = max(abs(uxqedge(:,q)));
%rms
uxrmstot(q) = sqrt(nanmean(uxqtot(:,q).^2));
uxrmscent(q) = sqrt(nanmean(uxqcent(:,q).^2));
uxrmsedge(q) = sqrt(nanmean(uxqedge(:,q).^2));

q

end









