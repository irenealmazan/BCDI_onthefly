% Matlabphasing_ver1.1.m 
%3D phase retrieval - Jesse Clark, LCN, UCL October-November 2010
%                       jesse.clark@ucl.ac.uk, jessenclark@gmail.com
clear
params.version='Matlab phasing version 1.1 - Nov 2010';
%%%%%%%%%%%%%%%
%% Select the data directory and the data files to load
%% The background files don't have to be supplied if this has already been done
%% Also select the size of the output and the binning and data threshold 
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

data_dir=dir;                   %assumes data is in the same dir as this file
                                %otherwise specify the data directory.
% at the moment it will create a save dir with a name that is 
% generated (see bottom of page) according to phasing params.  it will be
% saved in the directory that this phasing file is in.  this can be over-
% ridden by simply specifiying another save_dir.  
save_dir=dir;

% the data files and background files to load.  if no bg files are needed
% put [] in place of full_bg in bin_crop_center. i.e
% data=bin_crop_center(full_files,[],bin,min_data,aliens,nnc);

%load the crystals
name_of_this_file='Matlabphasing_ver1_1.m';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

load([dir,'crystal.mat'])

% create the sample
amp=zero_pad(zero_pad(crystal,250,250),256,256);
amp=round(amp/max(amp(:))*100)/100;
amp=imresize(amp,0.5,'nearest');
nn=size(amp);
sy=73*(nn(1)/256.0);  %largest dims of object
sx=80*(nn(2)/256.0);

px=.50;
py=.50;
[xx yy]=meshgrid(-(nn(1)/2-0.5):(nn(1)/2-0.5),-(nn(1)/2-0.5):(nn(1)/2-0.5));
xx=xx/max(abs(xx(:)));
yy=yy/max(abs(yy(:))); %norm so we can do it in terms of pi shifts

phase=exp(-i*pi*px*xx.^2/(sx/nn(2))^2-i*pi*yy.^2*py/(sy/nn(1))^2);

name=[num2str(px),'x-',num2str(py),'y-'];

% rsd=5;
% s = RandStream('mcg16807','Seed',rsd);
% RandStream.setDefaultStream(s)
% phase=exp(-i*convnfft(imresize(random('unif',0,2*pi,[20,20]),[nn(1),nn(2)],'nearest'),gauss_2D(11,11,1,1),'same'));

sample=amp.*phase;
params.hist=imhist(amp/max(amp(:)),100);
% creta data

data=abs(fftshift(fftn(fftshift(sample))));
nn=[size(data,2),size(data,1),size(data,3)];    %put the dimensions into xyz format
bin=[1,1];
comments='';
seq_let='';
files='';
back='';
aliens='';
min_data='';
save_data='';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Experimental Geometry
% used for displying the .rec files in vtk and for converting the .rec
% files to vtk

params.det_px=175.0e-6;              %detector pixel size (m)
params.lam=.13933;                   %wavelength (nm)
params.delta=32.300;                  %delta (degrees)
params.gam=10.80;                     %gamma (degrees)
params.arm=1.500;                    %camera distance (m)
params.dth=(1.2/61.0);              %dth angular range / number of images
params.dtilt=0;

params.spx=params.arm*params.lam/nn(1)/bin(1)/params.det_px;
params.spy=params.arm*params.lam/nn(2)/bin(2)/params.det_px;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER PARAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Algorithm selection and control      

start_guess='flat'; %start guess for iterative procedure. 'flat' uses the support
                    % 'random' uses a random object the size of the
                    % support.

                  %select one or two algorithms.  Choices are 'ER','SF','HIO','DM',
                  %'RAAR','HPR','ASR'.  Switching between algorithms is determined
                  %by trigger. to use one algorithm set, ALG2=ALG1
ALG1='ER' ;       %Algorithm 1.  must be in capitals
ALG2='HIO'  ;      %Algorithm 2.  must be in capitals

                    %trigger is an array of iteration numbers 
                    %that determine when to switch
                    %algorithms.  ie. trigger=[10,100,110,150] will do
                    %ALG1 until iteration 10 then switch to
                    %ALG2 until iteration 100 then it will
                    %switch back to ALG1 at 110 etc etc. until all iterations
                    %are complete.  There is no limit to how many times it
                    %can change.  The number of iterations is set by
                    %iterations though
                
beta=0.9;          %HIO,RAAR,HPR alg parameter
twin=0;             %set this to do another support (twin=1) that cuts out half the array
                    %it gets rid of the twin.  will only do it once at an
                    %iteration specified in the algorithm by twin_it.
twin_it=3;                    
%% 2. Support parameters, including shrinkwrap and phase constraint
% Parameters controlling the support.  To turn shrinkwrap OFF either set
% threshold=0 OR set the elements of shrink_trigger <0.

sx=0.35*nn(1);    %support sizes in x,y,z dimensions
sy=0.35*nn(2);    %i like to set them relative to the array dimensions   
sz=0.5*nn(3);    %but you can put in pixel values as well
 
threshold=0.00 ; %.1      %Shrinkwrap threshold (0<threshold<1).  threshold=0 will turn SW OFF.
sig=.9;                %Sigma of gauss kernal for smoothing in SW
sw_type='percent';      %Either 'gauss' or 'box' smoothing, with 'box' sig becomes
                      %the box size for smothing. [need to ceck this]
                              
shrink_trigger=sort([10:2:1000]);%[1,70]  ;        %the same as trigger but for controlling 
                                %shrinkwrap.  make all <0 to keep shrinkwrap 
                                %off or set threshold=0.
                                
phase_range=[-pi/2,0];          %phase range for phase constrained algorithms
phase_trigger=-[10];              %same rules as for trigger and shrink_trigger

%% 3. Partial coherence correction parameters                                 
pcdi=0;                          %do PCDI correction, 1 is ON, 0 is OFF

x_not_do=0             ;         %set = 1 to not do x
y_not_do=0           ;           %set = 1 to not do y
z_not_do=0 ;                     %set = 1 to not do z

angle_not_do=1          ;        %set = 1 to not do angle

roi=[16,16,16];                  %ROI pixel values to do pcdi min over
                                 %don't be stingey, 
                                 %more pixels=better=slower

pcdi_type='lucy'   ;         %type of pcdi correction. 
                                 % 'gauss' is the only one 
                                 %supported one at the moment.clea
                                
update_its=25;                    %every X iterations will update the kernal, saves time
start=10;                         %start after this many iterations
stop=10000;              %when to stop the updating of the kernal
pcdi_norm=1;                     %normalise before doing minimisation, 1=yes,0=no
pcdi_ER=0;                       %do normalisation on FT[ PiS[pn] ], that 
kernalxy=[0.5,0.5,.5];           %guess values for intial kernal, sigma [x,y,z].  not tat important
                                 %just don't make it too large.  ~.5 for
                                 %each is a good start guess
angle=0;                        %initial guess at angle
angle_start=0;                  %iteration

max_sigmas=[1.5,1.5,1.5]*3.0/bin(1);       %upper limit for the gauss kernal sigma size
                                %used in the minimisation.  only used when
                                %1) gaussian model is used AND 2) the
                                %minimisation procedure is constrained.
                                %else it has no effect.  lower limit is 000
min_sigmas=[0,0,0.1];
%% 4. Miscellaneous pramaters
params.no_zero=0;                %let zero measured intensity values float, 0=no,1=yes
                                 %will include statitically regularized
                                 %amplitude constraint in future versions
params.GPU=0;                    %will be used to trigger GPU acceleration
                                 %not used yet
params.charge=0;
params.silent=0;

%% 5. Genetic algorithm paramters
params.GA=0;                        %turns it on
params.PSO=0;

iterations= 10;   %Iterations to perform.
trigger=[5,560];%,160,201,360];%110,150];%,60,100];
params.generations=1;              %do this many generations  
params.population=1000;               %with this many individuals
params.silent=1;
params.pre_align=0;             %align to the 1st during iterative procedure to save time
                                %seems to effect th result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% other stuff.  no need to change anything.  makes a parameter structure
%% to save all the necessary info.  is saved at the end of the
%% reconstruction.
% used for storing all the parameters to output into a text file and some
% intial algorithm parameters
kern=gauss_3D(11,11,11,kernalxy(1),kernalxy(2),kernalxy(3));
params.comments=comments;
params.seq_let=seq_let;
params.files=files(:);
params.back=back;
params.aliens=aliens;
params.sz=size(data);
params.binning=bin;
params.min_data=min_data;
params.save_data=save_data;
params.iterations=iterations;
params.start_guess=start_guess;
params.ALG1=ALG1;
params.ALG2=ALG2;
params.trigger=trigger;

params.beta=beta;
params.size_init=[sx,sy,sz];
params.threshold=threshold;
params.sw_type=sw_type;
params.sigma=sig;
params.shrink_trigger=shrink_trigger;
params.phase_range=phase_range;
params.phase_trigger=phase_trigger;
params.twin=twin;
params.twin_it=twin_it;

params.pcdi=pcdi;
params.pcdi_type=pcdi_type;
params.update_its=update_its;
params.start_pcdi=start;
params.stop_pcdi=stop;
params.pcdi_norm=pcdi_norm;
params.pcdi_ER=pcdi_ER;
params.kernalxy=kernalxy;
params.angle=angle;
params.do=[x_not_do,y_not_do,z_not_do];
params.angle_not_do=angle_not_do;
params.angle_start=angle_start;
params.roi=roi;
params.max_sigmas=max_sigmas;
params.min_sigmas=min_sigmas;
params.orig=params.pcdi;
params.coh=kern;                 %guess at kerna
params.update=1;                 %tells the alg weather to update the coh 
                                 %func, no need to change this
ss=round([sy,sx,sz]);          %its row column major so y is before x    
if numel(size(data)) == 3
    support=zero_pad_ver2(ones(ss),nn(1),nn(2),nn(3) );    %create support 
else  support=zero_pad_ver2(ones([ss(1),ss(2)]),nn(1),nn(2) );end

if strcmp(start_guess,'flat'), pn=support;end%
if strcmp(start_guess,'random'), pn=support.*random('uniform',0,1,[nn(2),nn(1),nn(3)]);end%
%% Optimise the FFTW

fftw('planner', 'patient');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The phasing algorithm part.  itertive phasing handles the most part.
%% it calls modulus_projector and next_iterate to to the phasing.
%% modulus_projector contains all the pcdi stuff via determine_coh

%% GA
% seed
rsd=1;
s = RandStream('mcg16807','Seed',rsd);
RandStream.setDefaultStream(s)
% setup intial pn
npop=params.population;
ngens=params.generations;
iterations=params.iterations;
pn_g=zeros([params.sz,npop]);    %create empty array
supports=zeros([params.sz,npop]);
chi_ga=zeros([ngens,npop,iterations]);

L1_norm=zeros([ngens,npop]);

if ndims(data) == 3,
   for qq=1:npop
      pn_g(:,:,:,qq)=support.*random('uniform',0,1,[params.sz]); 
      supports(:,:,:,qq)=support;
   end
else
   for qq=1:npop
      pn_g(:,:,qq)=support.*random('uniform',0,1,[params.sz]); 
      supports(:,:,qq)=support;
   end
end


velocity=zeros(size(pn_g));
params.psi2=1.0;

pn_best=zeros(size(pn_g));
chi_best=zeros([npop,1]);
group_best=zeros(size(support));

%
for ww=1:ngens                %di ot over the generations, after each is where
                              %the 'breeding' takes place
    ref=0;
    counter=0;                          %counter to tell how many were updated
        
    disp(['GENERATION - ',num2str(ww)])
    
    for qq=1:npop             %iterate on each of the population

        if ndims(data) == 3,
            support=supports(:,:,:,qq);  %use the last support
            pn=pn_g(:,:,:,qq);end       %get the next iterate
        if ndims(data) == 2,
            support=supports(:,:,qq);
            pn=pn_g(:,:,qq);end

        %do the iterations
        [pnm support params_out chi DM_error coh ] = iterative_phasing(pn,data,support, params);

        if params.pre_align == 1
            if qq == 1,ref=pnm;end
            pnm=align_and_zero(ref,pnm); 
        end
        
        if ndims(data) == 3,
            supports(:,:,:,qq)=support;     %save new support
            pn_g(:,:,:,qq)=pnm.*support;end          %save new estimate
        if ndims(data) == 2,
            supports(:,:,qq)=support;
            pn_g(:,:,qq)=pnm.*support;end
        
        chi_ga(ww,qq,:)=chi;                % save chi
        L1_norm(ww,qq)=sum(abs(pnm(:)));
        
        
        if ww ~= 1,
            % check if the current is better, if then update best
            if chi_ga(ww,qq,end) < chi_best(qq),
                % update chi best
                chi_best(qq)=chi_ga(ww,qq,end);
                if ndims(data) == 2,pn_best(:,:,qq)=pn_g(:,:,qq);end
                if ndims(data) == 3,pn_best(:,:,:,qq)=pn_g(:,:,:,qq);end
                disp(['------UPDATING BEST [',num2str(qq),']------'])
                counter=counter+1;
            end
        end
        
       
        
    end

    if ww ~= 1
        if counter == 0,
            disp('NONE UPDATED')
            parama.psi2=params.psi2+1.0;
            disp(['INCREASING PSI2 - [',num2str(params.psi2),']'])
        end
    end
    
   
    
    % set the first iterate to the best since we don't have one yet
    if ww == 1,
        chi_best(:)=chi_ga(ww,:,end);
        pn_best=pn_g;
        disp('------SETTING FIRST ITERATE TO BEST------')
    end
    disp(squeeze(chi_best))
    % find the group best
    ind=find(chi_best == min(chi_best(:)));
    % chi best and pn_best will have the group best information and iterate
    if ndims(data) == 2,group_best=pn_best(:,:,ind);end
    if ndims(data) == 3,group_best=pn_best(:,:,:,ind);end
    disp(['------GROUP BEST IS [',num2str(ind),']------'])
 
    
%     
%     if ww >=9 & ww < 13
%         %disp('******SWITCHING TO PSO*******')
%         params.PSO=1;
%         params.GA=0; 
%     else
%         params.PSO=0;
%         params.GA=1;
%     end
    
    % PSO
    if params.PSO == 1
        if ww ~= ngens
            disp('------PSO UPDATING ITERATES------')
            [pn_g velocity]=PSO_iterates(pn_g,pn_best,velocity,group_best,chi_best,params);
        end
    end
    
    % GA
    if params.GA == 1
    %do the breeding
        if ww ~= ngens
            disp('------GA UPDATING ITERATES------')
            [ pn_g ] = GA_breed_iterates(pn_g,chi_ga(ww,:,end),params);%chi_ga(ww,:,end)
        end
       
        if threshold ~= 0
            supports=generate_new_supports(pn_g,params);
        end
    end
    
    
end

chi=chi_ga(end,:,end);
ind=find(chi == min(chi(:)));

pn=pn_g(:,:,ind);

%% Save the matlab reconstruction and output matlab files to be read
%% by mayavi.  outputs a phasingparams.py file and Matlab_to_vtk.py file

save([save_dir,name,'pn_g.mat'],'pn_g')
save([save_dir,name,'pn_best.mat'],'pn_best')

save([save_dir,name,'chi_ga.mat'],'chi_ga')

%save_matlab_phasing(pnm,support,data,coh,chi,params_out,save_dir  );
1
%% Display some things

%display_rec2D(pn,support,chi,DM_error);
%clear

%% END OF FILE %% END OF FILE 

1;