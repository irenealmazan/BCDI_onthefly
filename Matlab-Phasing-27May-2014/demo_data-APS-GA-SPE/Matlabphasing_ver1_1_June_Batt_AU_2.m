clc
clear all
params.version='Matlab phasing version 1.1 - Nov 2013';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% IN SITU June DATA PARAMS, for doing single scans

% snums=[196 197 199 200 203 205 207 209 210 212];
% 
% snums = 1072;

%snums=[705 717 735 750 761 772 780 781 791 801 840 850 854];

for tt = 1
	
for tq = 1


%% Determining the current directory
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data preperation
params.data_dir={'/Users/a4894z/Dropbox/Matlab_phasing_AUG13/GA/'};            %assumes data is in the same dir as this file
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

%params.files={[int2str(snums(tt)) '_a.mat']};
params.files={'196.tif'};
%params.files={'705_a.mat'};%,'1942.tif','1943.tif'};%,'1924.tif','1925.tif'};

params.back={};%if no bg required, leave empty {}.  don't use {''}.

params.itnum = 1;               %change this for each recon of same data set
 
params.seq_let=['r' num2str(tq)]; %sequence letter, second part of save name

%params.comments2 = ['Novdef' num2str(snums(tt))];
params.comments2 = ['test'];

params.comments='medium support, coord transform';      %write comments here and this will be saved in the params file

params.binning=[1,1];           %binning for x and y

params.skipping=0 ;             %will bin (=0), will skip points (=1)

params.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=1;            %min data threshold.  below this is removed.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.
                                
params.schot_th=0;           %secondary threshold applied AFTER binning
params.subtract_dc=0;           %leave this as 0 (unless you know what you are doing)
params.no_center=0;             %no centering of data -> no_center=1
params.no_fft_pad=0;            %no padding for fft -> no_fft_pad=1
params.pad_ptych=0;             %leave this as 0 (unless you know what you are doing)
params.no_hist=1;               %plot the histograms of the data (=0)
params.bg_mult=1;               %multiplication to apply to the bg file (if exp time is different)
params.save_data='YES';


params.data_only='NO';        
                                %will save or exit after data prep. ='NO' no save and no exit,         
                                %'YES' will save and continue, 'YES-EXIT'
                                %will save and then stop the script

params.nnc=[0,0,0,0,5,5];
                                % initial cropping of data before binning.
                                %eg. nnc=[0,0,-10,-10,5,5] will do nothing to x,
                                %will crop 10 pixels off each end in y and will pad
                                %5 pixels to each end in z.  set nnc=0 to do
                                %nothing.
params.do_2D=0;                 %will take the central slice and do 2D phase retrieval if =1                                
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% No inputs required here.  Doing the data preperation.                     
[data params] = load_MP_data(params);   %load the data


%% check for saturated pixels
  [p1 p2 p3] = ind2sub(size(data),find(data>=1000));
   %condition of max being greater than 150% of n.n.
   %length(p1)
   for qq=1:length(p1)
   disp('checking for saturated pixels')
   nn1 = data(p1(qq)+1,p2(qq),p3(qq));
   nn2 = data(p1(qq)-1,p2(qq),p3(qq));
   nn3 = data(p1(qq),p2(qq)+1,p3(qq));
   nn4 = data(p1(qq),p2(qq)-1,p3(qq));
   mval = data(p1(qq),p2(qq),p3(qq));
   
   if mval>500 && mval>=1.5*nn1 && mval>=1.5*nn2 && mval>=1.5*nn3 && mval>=1.5*nn4
  %     disp('saturated pixel found with value')
  %     mval
  %     disp('changing value to')
       newval = (1/4)*(nn1+nn2+nn3+nn4);
 %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
 %      pause
       data(p1(qq),p2(qq),p3(qq))=newval;
 %      scannum
 %      close(gcf)
   end
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Experimental Geometry
% used for displying the reconstructions

params.det_px=55e-6;               %detector pixel size (m)
params.lam=.1379;                   %wavelength (nm)
params.delta=17.0135;                  %delta (degrees)
params.gam=0;                     %gamma (degrees)
params.arm=1.823175;                      %camera distance (m)
params.dth=.01;                      %angular step size
params.dtilt=0;


params.spx=params.arm*params.lam/params.nn(1)/params.binning(1)/params.det_px;   %sample pixel size
params.spy=params.arm*params.lam/params.nn(2)/params.binning(2)/params.det_px;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Algorithm selection and control      
params.iterations= 2150;             % Iterations to perform.

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
params.ALG1='HIO' ;                  %Algorithm 1.  must be in capitals
params.ALG2='ER' ;                 %Algorithm 2.  must be in capitals

%params.trigger=[10,40,50,80,90,120,130,160,170,200,210,240,250,280,290,320,330,360,370,400,410,440,450,480,490,520,...
%	530,560,570,600,610,640,650,670,680,690,710,720,730,750,760,770,780,800];

params.trigger=sort([(1:90:params.iterations),(1:90:params.iterations)+10,params.iterations-10]);

%params.trigger(1) = 5;


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



% sx=0.2*params.nn(1);                  %x,y,z support size (pixels)
% sy=0.2*params.nn(2);                  %i like to set them relative to the array dimensions   
% sz=0.35*params.nn(3);                  %but you can put in pixel values as well

sx=80;sy=80;sz=50;

 
params.threshold=0.1;                 %Shrinkwrap threshold (0<threshold<1).  threshold=0 will turn SW OFF.
params.sigma=1;                         %Sigma of gauss kernal for smoothing in SW
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
                                       
                              
params.shrink_trigger=generate_trigger(5,4,params.iterations);
                                       %the same as trigger but for controlling 
                                       %shrinkwrap.  make all <0 to keep shrinkwrap 
                                       %off or set threshold=0.
                                
params.phase_range=[-pi/2,pi/2];       %Pixels with a phase outside this range will have there
                                       %amplitude->0. Done within the
                                       %initial support region or
                                       %shrinkwrap support region
                                       
params.phase_trigger=[1,10];         %same rules as for trigger and shrink_trigger
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
                                 
params.update_its=25;               %every X iterations will update the kernal, saves time
params.start_pcdi=150;               %start doing PCDI correction after this many iterations
params.stop_pcdi=params.iterations; %when to stop the updating of the kernal
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

params.GA_metric='sharpness';           %how to decide what is the 'best'
                                        %us used for the breeding only
                                        %'chi' - usual chi squared (min
                                        %chi)
                                        %'sharpness' - {|p|^4 (min
                                        %sharpest)
                                        %'area' -  the area of the support (largest) 
                                        %'TV' - min of the TV norm
params.GA_metric_return='sharpness';    % how to decide what to return at the end
                                            
                                        
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
params.population=5;                                      
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
params.iterate_avg_start=1000;    %when to start averaging iterates
params.iterate_avg_its=50;        %number of iterates to average
params.iterate_avg_inter=100;      %interval between adding averages

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


%%decide if you want to use a box support or the previous amplitudes

%box support
[pn support] = make_intial_support(data,ss,params);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The phasing algorithm part.  itertive phasing handles the most part.
%% it calls modulus_projector and next_iterate to to the phasing.
%% modulus_projector contains all the pcdi stuff via determine_coh
[pnm support params_out chi DM_error coh] = iterative_phasing_general(pn,data,support, params);

%% Save the matlab reconstruction and output matlab files to be read
%% by mayavi.  outputs a phasingparams.py file and Matlab_to_vtk.py file

params_out.newxyz=[256,256,128];

save_matlab_phasing(pnm,support,data,coh,chi,params_out,params.save_dir);

%% Display some things

display_rec3D(pnm,support,chi,DM_error);
%clear
close all;

tt

%% END OF FILE %% END OF FILE 

end

beep

end


%%
%Need to load, check for twins, and average the resulting reconstructions
clear all;
numscans=5;
%some scans are 64 and some are 48 so be aware
%iterates = zeros(256,256,48,numscans);
%cohers = zeros(size(iterates));
%chis = zeros(1,numscans);
linup = 1; %try to shift to match



for snum=[207]

%do the loading into a 4D matrix
for i=1:numscans

datdir = ['/Users/andrewulvestad/Desktop/J' int2str(snum) 'r' int2str(i) '/'];

%fname = ['J' int2str(snum) 'r' int2str(i) '-LAB.mat']; %phase and amp

%fname2 = ['J' int2str(snum) 'r' int2str(i) '-LAB-S.mat']; %support

%fname3 = ['J' int2str(snum) 'r' int2str(i) '-LAB-C.mat']; %coherence

fname4 = ['J' int2str(snum) 'r' int2str(i) '-ERROR.mat']; %error, just take final value

fname5 = ['J' int2str(snum) 'r' int2str(i) '-AMP.mat'];

fname6 = ['J' int2str(snum) 'r' int2str(i) '-PH.mat'];


%anew = load([datdir fname]); %load amp/phase
%cnew= load([datdir fname3]); %load coh
%snew = load([datdir fname2]);
enew = load([datdir fname4]);
amps = load([datdir fname5]);
phs = load([datdir fname6]);

%anew = anew.array;
% cnew = cnew.array;
% snew = snew.array;
% enew = enew.array;
amps = amps.array;
phs = phs.array;

[xlim ylim zlim] = size(amps);


%amps(amps>0.3)=1;

%figure;

%imagesc(amps(:,:,zlim/2)); zoom(4);

%use the nontransformed ones
iterates(:,:,:,i) = amps.*exp(1i.*phs);

chis(i) = enew.chi(end);

%cohers(:,:,:,i) = cnew.array;

%supports(:,:,:,i) = snew.array;


end

% [ylim xlim zlim] = size(anew.array);
% 
% %only keep certain ones
% excld = [10 8 6 1];
% 
%  iterates(:,:,:,excld)=[]; 
%  %cohers(:,:,:,excld)=[];
%  supports(:,:,:,excld)=[];
%  chis(:,excld)=[];
 
 [x1 y1 z1 q1] = size(iterates)

%now do the averaging
ref = iterates(:,:,:,1);

for qq=2:q1   
        
   
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



asum = mean(iterates,4); %average amp and phase at same time
%csum = mean(cohers,4);
%ssum = mean(supports,4);
chit = mean(chis)

end

%%
asum = mean(iterates,4);
smamp = smooth3(abs(asum),'box',[7 7 3]);
figure; isosurface(smamp,0.05);
%%
nm=5;
for i=z1/2-10:z1/2+20
%	figure;
%	imagesc(angle(iterates(:,:,i,nm)).*(supports(:,:,i,nm)>=1)); caxis([-pi pi]); zoom(2); axis off;
     %figure;
 	
     %imagesc(angle(iterates(:,:,i,nm)).*(ssum(:,:,i)>=1));
		 
		 %phasesp = angle(iterates(:,:,:,nm)).*(smamp>=0.05);
         phasesp = angle(asum).*(smamp>=0.05);
    
    imagesc(phasesp(:,:,i));
    
		
    caxis([-pi pi]); zoom(2); %axis off;
    %title(['phases for scan ' num2str(nm) ' at height ' num2str(i)])
pause
end

%%
%get the strain first
phasesp = angle(asum).*(abs(asum)>=0.2);
phasesp2 = phasesp;
phasesp2(phasesp2==0)=nan;

[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom

[ux2 uy uz] = gradient(phasesp2,160,160,470);

%coordinate invariant metric
Estrain = nansum(nansum(nansum(ux2.^2)))*((16e-9)*(16e-9)*(47e-9))


%put all these bad boys into a structure

%TOT.asum = asum;
%TOT.csum = csum;
%TOT.ssum = ssum;
%TOT.chit = chit

v43 = genvarname(['Jpart' int2str(snum)]);
eval([v43 '= asum;']);

%make the images for the excel sheet
% subplot(2,1,1)
% imagesc(abs(asum(:,:,zlim/2))); zoom(2)
% subplot(2,2,2)
% imagesc(angle(asum(:,:,zlim/2))); zoom(2); caxis([-pi pi]);
figure;
imagesc(ux(:,:,zlim/2).*(abs(asum(:,:,zlim/2))>=0.2)); caxis([-1e-3 1e-3]) ; zoom(2)

figure;
isosurface(abs(asum),0.2)

%save the average in a mat file

save(['Jpart' int2str(snum)],v43); %saves it as 40A.mat, x40A is asum

clear TOT asum csum ssum chit iterates cohers supports chis;


%%
%assuming you've loaded them all in and they're all in structures
clear iterates;
i=1;

for snum=[45 46]

v43 = genvarname(['Jpart' int2str(snum)]);
eval(['TOT=' v43;]);

[ylim xlim zlim] = size(TOT.asum);

%enlarge to 64 if you're not there
if zlim~=64
    TOT.asum=padarray(TOT.asum,[0 0 (64-zlim)/2]); 
end



%compute the strain and check that out
phasesp = angle(TOT.asum).*(abs(TOT.asum)>=0.2); %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp/(2*pi*sqrt(3)/8.14),160,160,350); %use spacing of 20nm or 200 angstrom
ux(isnan(ux)==1) = 0;

%figure generation, decide what you want
% cdata = smooth3(ux,'box',[3 3 1]);
% %subplot(5,2,i); 
%   figure;
%   imagesc(cdata(:,:,zlim/2));
%   %caxis([-20e-5 20e-5]); 
%   colorbar;
%    zoom(4);
 
% figure;
% [x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
% cdata = smooth3(ux,'box',[3 3 1]);
% p = patch(isosurface(x,y,z,abs(TOT.asum),0.2));
% isonormals(x,y,z,abs(TOT.asum),p);
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

%look at the phases but cut up
% figure;
% imagesc(angle(TOT.asum(:,:,64/2)).*(abs(TOT.asum(:,:,64/2))>=0.35))




% figure;
% imagesc(ux); zoom(5)

%figure;
%isosurface(abs(TOT.asum),0.25)

figure; imagesc(abs(TOT.asum(:,:,64/2)));

iterates(:,:,:,i) = TOT.asum;

i=i+1;

end

[x1 y1 z1 q1] = size(iterates)
ref = iterates(:,:,:,1);
linup=1;

for qq=2:q1     
        
   
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

%check what they look like after lineup and cnj_refl

for qq=1:q1
figure;
imagesc(abs(iterates(:,:,zlim/2,qq))); zoom(4);
end


%now we can average together and get the amplitude
asum = mean(abs(iterates),4);

for i=1:64
    imagesc(asum(:,:,i))
    colorbar;
    pause
end

%and look at phase info 
for itnum=1:2
%cut with individual
%phasesp = angle(iterates(:,:,:,itnum)).*(abs(iterates(:,:,:,itnum))>=0.4); %this is now displacement in angstrom
%cut w average
phasesp = angle(iterates(:,:,:,itnum)).*(asum>=0.5);

phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,200,200,350); %use spacing of 20nm or 200 angstrom


% figure;
% imagesc(angle(iterates(:,:,zlim/2,itnum)))


%just look at slices
 figure;
cdata2 = smooth3(ux,'box',[3 3 1]);
 imagesc(cdata2(:,:,zlim/2)); zoom(4)



%strain projected onto an isosurface
% figure;
% ux(isnan(ux)==1) = 0;
% cdata = smooth3(ux,'box',[3 3 1]);
% 
% [ylim xlim zlim] = size(asum);
% [x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
% p = patch(isosurface(x,y,z,asum,0.5));
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
caxis([-2e-4 2e-4]); 
% %caxis([-1 1]);
% set(gca,'fontsize',16)

end



%% if you're doing GA then you can use this. just load in the Lab.mat file

clear all;

i=1;

for snum=1902

%do the loading into a 4D matrix
datdir = ['/Users/andrewulvestad/Desktop/recon/J' int2str(snum) 'GAr' int2str(1) '/'];

fname = ['J' int2str(snum) 'GAr' int2str(1) '-LAB.mat']; %phase and amp

fname2 = ['J' int2str(snum) 'GAr' int2str(1) '-LAB-S.mat']; %support

fname3 = ['J' int2str(snum) 'GAr' int2str(1) '-LAB-C.mat']; %coherence

fname4 = ['J' int2str(snum) 'GAr' int2str(1) '-ERROR.mat']; %error, just take final value

fname5 = ['J' int2str(snum) 'GAr' int2str(1) '-AMP.mat'];

fname6 = ['J' int2str(snum) 'GAr' int2str(1) '-PH.mat'];


anew = load([datdir fname]); %load amp/phase
cnew= load([datdir fname3]); %load coh
snew = load([datdir fname2]);
enew = load([datdir fname4]);
amps = load([datdir fname5]);
phs = load([datdir fname6]);

anew = anew.array;
% cnew = cnew.array;
% snew = snew.array;
% enew = enew.array;
amps = amps.array;
phs = phs.array;

[xlim ylim zlim] = size(anew);

%figure;

%imagesc(amps(:,:,zlim/2)); zoom(4);

%amps(amps>0.3)=1;

%figure;

%imagesc(amps(:,:,zlim/2)); zoom(4);

%use the nontransformed ones
iterates(:,:,:,i) = amps.*exp(1i.*phs);

asum = amps.*exp(1i.*phs);

[xlim ylim zlim] = size(asum);

chis(i) = enew.chi

i=i+1;

phasesp = angle(asum).*(abs(asum)>=0.2);
phasesp2 = phasesp;
phasesp2(phasesp2==0)=nan;

[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom

[ux2 uy uz] = gradient(phasesp2,160,160,470);

%coordinate invariant metric
Estrain = nansum(nansum(nansum(ux2.^2)))*((16e-9)*(16e-9)*(47e-9))

%make the images for the excel sheet
subplot(2,2,1)

imagesc(abs(asum(:,:,zlim/2))); zoom(2)
subplot(2,2,2)
imagesc(angle(asum(:,:,zlim/2))); zoom(2); caxis([-pi pi]);
cdata2 = smooth3(ux,'box',[3 3 1]);
subplot(2,2,3)
imagesc(cdata2(:,:,zlim/2)); caxis([-1e-3 1e-3]) ; zoom(2)
subplot(2,2,4)
isosurface(abs(asum),0.2)
title(int2str(snum))


end

 
[x1 y1 z1 q1] = size(iterates)



for qq=[2 4]     
iterates(:,:,:,qq)=conj_reflect(iterates(:,:,:,qq));
end


%now do the averaging
ref = iterates(:,:,:,1);



for qq=2:q1   
        
    [h k l] = register_3d_reconstruction(abs(ref),abs(iterates(:,:,:,qq)));
    next_c = sub_pixel_shift(iterates(:,:,:,qq),h,k,l);
    iterates(:,:,:,qq) = next_c;

    
end


%what do they look like after all this?
for qq=1:q1
subplot(2,2,qq)
imagesc(abs(iterates(:,:,zlim/2,qq)));
end

for qq=1:q1
subplot(2,2,qq)
imagesc(angle(iterates(:,:,zlim/2,qq))); caxis([-pi pi])
end

%can check phase but no gurantee they get the same global phase offset so
%check strain instead
for qq=1:q1
figure;
imagesc(angle(iterates(:,:,zlim/2,qq)).*(abs(iterates(:,:,zlim/2,qq))>=1));
caxis([-pi pi])
end


for qq=1:q1

phasesp = angle(iterates(:,:,:,qq)).*(abs(iterates(:,:,:,qq))>=0.2);
phasesp2 = phasesp;
phasesp2(phasesp2==0)=nan;

[ux uy uz] = gradient(phasesp,200,200,350); %use spacing of 20nm or 200 angstrom

[ux2 uy uz] = gradient(phasesp2,200,200,350);

%coordinate invariant metric
Estrain(qq) = nansum(nansum(nansum(ux2.^2)))*((16e-9)*(16e-9)*(47e-9));

%just look at slices
 figure;
cdata2 = smooth3(ux,'box',[3 3 1]);
 imagesc(cdata2(:,:,zlim/2)); caxis([-1e-3 1e-3]) 
 
end

%% load up one scan for each reconstruction

tv=1;

for snum=[196 197 199 200 203 205 207 209 210 212]

%do the loading into a 4D matrix
for i=1

datdir = ['/Users/andrewulvestad/Desktop/recon/JC' int2str(snum) 'r' int2str(i) '/'];

%fname = ['J' int2str(snum) 'r' int2str(i) '-LAB.mat']; %phase and amp

fname2 = ['JC' int2str(snum) 'r' int2str(i) '-SUP.mat']; %support

%fname3 = ['J' int2str(snum) 'r' int2str(i) '-LAB-C.mat']; %coherence

fname4 = ['JC' int2str(snum) 'r' int2str(i) '-ERROR.mat']; %error, just take final value

fname5 = ['JC' int2str(snum) 'r' int2str(i) '-AMP.mat'];

fname6 = ['JC' int2str(snum) 'r' int2str(i) '-PH.mat'];


%anew = load([datdir fname]); %load amp/phase
%cnew= load([datdir fname3]); %load coh
snew = load([datdir fname2]);
enew = load([datdir fname4]);
amps = load([datdir fname5]);
phs = load([datdir fname6]);

%anew = anew.array;
% cnew = cnew.array;
snew = snew.array;
% enew = enew.array;
amps = amps.array;
phs = phs.array;

[xlim ylim zlim] = size(amps);

%figure;

%imagesc(phs(:,:,zlim/2)); zoom(4);

%amps(amps>0.3)=1;

%figure;

%imagesc(amps(:,:,zlim/2)); zoom(4);

%use the nontransformed ones
iterates(:,:,:,tv) = amps.*exp(1i.*phs);

chis(i) = enew.chi(end);

%cohers(:,:,:,i) = cnew.array;

supports(:,:,:,tv) = snew;



end
tv=tv+1;

end

%% line them up

ref = iterates(:,:,:,1);

for qq=2:tv-1  
        
   
    cnj_rf=is_conj_ref(ref,iterates(:,:,:,qq))
            
    if cnj_rf ==1
        iterates(:,:,:,qq)=conj_reflect(iterates(:,:,:,qq));
    end 
    
    %try to line them up using essentially cross correlation
    
    [h k l] = register_3d_reconstruction(abs(ref),abs(iterates(:,:,:,qq)));
    next_c = sub_pixel_shift(iterates(:,:,:,qq),h,k,l);
    iterates(:,:,:,qq) = next_c;

    
end

asum2 = mean(abs(iterates),4);
ssum = mean(supports,4);
pval=0.1;

supavg = mean(squeeze(mean(ssum,1)),1);

p3 = find(supavg==max(supavg));

%check what they look like after that
for qq=1:tv-1
    figure; imagesc(abs(iterates(:,:,p3,qq)))
    pause
end



%% lattice constant vector (taken from excel file)

%for june exp 4 post cycling
latc = [8.173212187
8.179662459
8.175844893
8.169025791
8.114736231
8.098563875
8.143313028
8.174397507
8.193170026
8.200700122];

voltc = [3.57
4.08
4.42
4.68
4.74
4.68
4.6
4.42
3.95
3.66];



figure; isosurface(abs(asum2),pval); 
figure; imagesc(abs(asum2(:,:,p3)))

%support stuff
figure; imagesc(ssum(:,:,p3))


%try just replacing the amplitudes with the support and then cnj refl




%% strain computation
for qr=1:tv-1
%redefine to each one of the iterates
anew = iterates(:,:,:,qr);
snew = smooth3(supports(:,:,:,qr),'box',[5 5 1]);

%phasesp = angle(anew).*(snew>=.9);

%phasesp = angle(anew).*(abs(anew)>=0.1);

phasesp = angle(anew);

 phasesp2 = phasesp;
 phasesp2(phasesp2==0)=nan;

[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom

 [ux2 uy uz] = gradient(phasesp2,160,160,470);

%need to take out the maximal strain values and replace by nearest
%neighbors

%[p1 p2 p3] = ind2sub(size(ux2),find(abs(ux2)>=.01));
   %condition of max being greater than 150% of n.n.
   %length(p1)
%    for qq=1:length(p1)
%    
%    nn1 = ux2(p1(qq)+1,p2(qq),p3(qq));
%    nn2 = ux2(p1(qq)-1,p2(qq),p3(qq));
%    nn3 = ux2(p1(qq),p2(qq)+1,p3(qq));
%    nn4 = ux2(p1(qq),p2(qq)-1,p3(qq));
%    mval = ux2(p1(qq),p2(qq),p3(qq));
%    
%    if  mval>=1.1*nn1 && mval>=1.1*nn2 && mval>=1.1*nn3 && mval>=1.1*nn4
%        disp('saturated pixel found with value')
%        mval
%        disp('changing value to')
%        newval = (1/4)*(nn1+nn2+nn3+nn4)
%  %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
%  %      pause
%        ux2(p1(qq),p2(qq),p3(qq))=newval;
%  %      scannum
%  %      close(gcf)
%  
% 	 elseif isnan(nn1)==1 | isnan(nn2)==1 | isnan(nn3)==1 | isnan(nn4)==1
% 		 newval = nan;
% 		 disp('saturated pixel found with value')
%        mval
%        disp('changing value to nan')
% 		 ux2(p1(qq),p2(qq),p3(qq))=newval;
% 		
% 	 end
% 	 
% 	
% 	 end
	 

%coordinate invariant metric
Estrain(qr) = nansum(nansum(nansum(ux2.^2)))*((16e-9)*(16e-9)*(47e-9));

%find the max of the strain within a small region
% mastr = max(max(max(ux.*(abs(anew>=0.5)))))
% mistr = min(min(min(ux.*(abs(anew>=0.5)))))

 %subplot(5,2,qr)
%   figure;
%   imagesc(ux(:,:,p3)); %caxis([mistr mastr]) ; 
%   caxis([-5e-3 5e-3]);
%   zoom(2); axis off; %colorbar;
%   
  figure;
  imagesc(phasesp(:,:,p3)); %caxis([mistr mastr]) ; 
  caxis([-pi pi]);
  zoom(2); axis off;

%subplot(3,3,qq)
%imagesc(abs(asum(:,:,s/2)))

%explore some of the individual scans amp/phase/strain
%figure;
%isosurface(ux2,0.5e-3)


strmap(:,:,:,qr) = ux; %finally store the strains

clear anew ux;

end
%% plots for strain and lat const and voltage

figure; plot(Estrain,'.-')
figure; plot(latc,'.-')
figure; plot(voltc,'.-')

%% investigate individual scans

figure;
nm=1;
for i=48-16:48+16
	imagesc(strmap(:,:,i,nm).*(ssum(:,:,i)>=0.95)); caxis([-1e-3 1e-3])
	pause
end

for i=48-16:48+16
	imagesc(ssum(:,:,i)); colorbar;
	pause
end



%% multiple contour slice method
for qr=1:10 %[3 5 6]
	
	ux2 = strmap(:,:,:,qr).*(ssum>=0.9);
	
	figure;
%this is a really great viz tool
%cdata = smooth3(ux2,'box',[3 3 1]);%.*(asum2>=0.3);
p3=96/2;
h = contourslice(ux2, [], [], [p3-3 p3-2 p3-1 p3 p3+1 p3+2 p3+3],256);
camva(24); camproj perspective;
campos([-3, -15, 5]);
set(gcf, 'Color', [.3, .3, .3], 'renderer', 'zbuffer');
set(gca, 'Color', 'black', 'XColor', 'black', ...
'YColor', 'black', 'ZColor', 'black');
box off;
view(-38,16)
zoom(3.5);
axis off;
set(gca,'color','k')
set(h, 'LineWidth', 6);
	
	
end


%% isosurf projection method
for qr=1:10
	
ux2 = strmap(:,:,:,qr);
	
figure;
%cdata = smooth3(ux2,'box',[3 3 1]).*(asum2>=0.3);
[x y z] = meshgrid(1:256,1:256,1:96);
cdata = smooth3(ux2,'box',[3 3 1]);
p = patch(isosurface(x,y,z,asum2,0.25));

p = patch(isosurface(x,y,z,abs(iterates(:,:,:,qr)),0.25));

isonormals(x,y,z,asum2,p);
isocolors(x,y,z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none')
%daspect([1 1 1]);
axis tight
camlight(230,250); 
lighting GOURAUD; 
%title('4.7 1st Charge ','FontSize',16)
colorbar;
%set(gca, 'ydir','reverse')
caxis([-1e-3 1e-3])
view(3); 
axis off;

	
end
	
	
	


