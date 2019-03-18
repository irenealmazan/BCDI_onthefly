function [pnm support params_out chi DM_error coh ] = iterative_phasing(pn,data,support, params)
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
% does the iterative phasing loop
% params file passes information baout when to change
%algorithms, support, pcdi correction etc etc
%returns 
tic
nn=params.sz;       %size of the data, sz=size(data)

if ndims(pn) == 2,params.do_2D = 1;end         %tell it it is doing a 2d rec (used for saving)

if params.twin == 1
    if ndims(pn) == 3
        support0=zeros([nn(1),nn(2),nn(3)]);            %support0
        support0(1:nn(1)/2,1:nn(2)/2,:)=1;
    end
    if ndims(pn) == 2
       support0=zeros([nn(1),nn(2)]);            %support0
       support0(1:nn(1)/2,1:nn(2)/2)=1;  
    end
        
else support0=1;end

support_orig=support;                   %keep original support for phase constraint stuff

params=set_param_defaults(params);

disp_off=params.silent;   %this is to turn the screen output off during the iterations
                          %params.silent=1 -> display off
                         
sw=0;
ph_cs=0;

chi=0;
DM_error=0;
pnm=0;
flip=1;                 %1 is alg 1
shrink_flip=-1;          %-1 is no shrink
phase_flip=-1;           %-1 is no phase constraint

gauss_f_flip=-1;        %-1 is no gauss filtering

try
    angle_orig=params.angle_not_do;
catch
    angle_orig=0;
end

if params.GPU == 1
    try
        support0=gdouble(support0);
    end
    support_orig=gdouble(support_orig);
    data=gdouble(data);
    support=gdouble(support);
    pn=gdouble(pn);
end

%check to see if low res to high res is activated and set some defualts
params=check_lres(params);

%%
%work if we want to save iterates during for GA type averaging
if params.iterate_avg_save == 1, 
    sharpness=[];
    pn_g=zeros([size(data),numel(params.iterate_avg_trigger)]);
end

pnm_avg=0;                       %will be used to store the average
avg_counter=0;
params.iterate_avg_enabled=0;    %marker for averaging, will append save name if =1

%% Iterative stuff
for qq = 1:params.iterations
    
    params.itno=qq;         %store current iteration number
    
    flip=calc_state(params.trigger,qq).*flip;       %check if it time to switch algs
    
    pn0=pn;
    
    if params.simultaneous_coh_det ~=1              %when doing the simultaneous recs, don't want to update coh function here
        if mod(qq,params.update_its) == 0, params.update=1;else params.update=0;end
        if qq >= params.start_pcdi, params.pcdi=params.orig; else params.pcdi=0;end 
        if qq >= params.stop_pcdi,params.update=0;end
    else params.update=0;end
        
    if qq < params.angle_start,params.angle_not_do = 0;else params.angle_not_do=angle_orig;end
    
    gauss_f_flip=calc_state(params.gauss_f_trigger,qq).*gauss_f_flip;
    params.gauss_f_flip=gauss_f_flip; %used for mod proj to turn on or off
    
    if disp_off ==0;disp('----------------------------------------------');end
    
    %check if doing a low res rec
    [ params ] = check_lres_it(params,flip,qq); %determine if ALG1 or ALG2 are low res
    sigxy=params.nonGA_lres_det(qq);            %get the detetctor mask width (std)
    
    if flip == 1    
        ALG=params.ALG1;
        switch params.ALG1
            case {'DM','DMr','DMlr','DMrlr'}
                [pnm error params]=modulus_projector(2*pn.*support-pn,data.*nonGA_dmask(params,sigxy),params,support);    
            otherwise
                [pnm error params]=modulus_projector(pn,data.*nonGA_dmask(params,sigxy),params,support);
        end
        %if strcmp(params.ALG1,'DM') || strcmp(params.ALG1,'DMr')
        %    [pnm error params]=modulus_projector(2*pn.*support-pn,data,params,support);
        %else 
        %    [pnm error params]=modulus_projector(pn,data,params,support);
        %end
    else
        ALG=params.ALG2;
        
        switch params.ALG2
            case {'DM','DMr','DMlr','DMrlr'}
                [pnm error params]=modulus_projector(2*pn.*support-pn,data.*nonGA_dmask(params,sigxy),params,support);    
            otherwise
                [pnm error params]=modulus_projector(pn,data.*nonGA_dmask(params,sigxy),params,support);
        end
        %if strcmp(params.ALG2,'DM') || strcmp(params.ALG2,'DMr')
        %    [pnm error params]=modulus_projector(2*pn.*support-pn,data,params,support);
        %else 
        %    [pnm error params]=modulus_projector(pn,data,params,support);
        %end
    end
    
    chi(qq)=error;
    
    shrink_flip=calc_state(params.shrink_trigger,qq).*shrink_flip;
    if params.threshold == 0,shrink_flip =0;end                     %override trigger if threshold=0
    phase_flip=calc_state(params.phase_trigger,qq).*phase_flip;
    
    if disp_off == 0,
        disp(['iteration number: (',num2str(qq),'/',num2str(params.iterations),')','[',ALG,']','  error =',num2str(error)])
    end
        
    if shrink_flip == 1   %doing it off pn worked better
        
        if strcmp(params.sw_type,'percent-auto') == 1
           
           disp('-------------------------------------------------------') 
           disp('| Calculating % area for support from autocorrelation....') 
           params.sw_type='gauss_percent';        %set back to percent
           
           auto_c=abs(fftshift(ifftn(fftshift(data.^2))));
           temp=shrink_wrap(abs(auto_c).^.5,.1,.2);  %.^.5,.05,.1
           params.threshold=sum(temp(:))/numel(data);
           
           if ndims(data) == 3,params.threshold=params.threshold/8.0;end
           if ndims(data) == 2,params.threshold=params.threshold/4.0;end
           
           clear temp
           clear auto_c
           disp(['| Using % of [',num2str(params.threshold*100.0),']'])
           disp('-------------------------------------------------------')
           
        end
       
        switch params.sw_type
            
            case 'std'
                params.threshold=std(abs(pn(:)))*params.sigma;
                support=shrink_wrap(abs(pn),params.threshold,.01,'gauss');
            
            case 'gauss-minarea'
                support=shrink_wrap(abs(pn),params.threshold,params.sigma,'gauss');
                if sum(support(:)) < params.sw_min_area*numel(support(:))
                    support=shrink_wrap(abs(pn),params.sw_min_area,params.sigma,'gauss_percent');
                    disp(['Minimum support area not met, setting to minimum of ',num2str(params.sw_min_area)])
                end
                
            case 'gauss-maxarea'
                support=shrink_wrap(abs(pn),params.threshold,params.sigma,'gauss');
                if sum(support(:)) > params.sw_max_area*numel(support(:))
                    support=shrink_wrap(abs(pn),params.sw_max_area,params.sigma,'gauss_percent');
                    disp(['Support area too large, setting to maximum of ',num2str(params.sw_max_area)])
                end
                
            case 'gauss-const'
                
            case 'gauss-pnm'    
                support=shrink_wrap(abs(pnm),params.threshold,params.sigma,'gauss');
            case 'gauss-square'
                support=shrink_wrap(abs(pn).^2,params.threshold,params.sigma,'gauss');
            otherwise
                support=shrink_wrap(abs(pn),params.threshold,params.sigma,params.sw_type);
               
        end
            
%         if strcmp(params.ALG1,'DMr') %| strcmp(params.ALG1,'HIOs')
%             support=shrink_wrap(abs(pnm),params.threshold,params.sigma,params.sw_type);
%         else
%             support=shrink_wrap(abs(pn),params.threshold,params.sigma,params.sw_type);
%         end
        
        if disp_off == 0,disp(['Shrinkwrap support, sig [',num2str(params.sigma),'] th [',num2str(params.threshold),']']),end

        sw=1;       %used to tell if sw was used at all
    else
        if disp_off == 0,disp('Fixed support'),end
    end
    
    if phase_flip == 1
        ph_support = phase_constraint(pnm,params.phase_range);
        
        if shrink_flip == 0
            support=support_orig.*ph_support;end
        
        if shrink_flip == 1
            support=support.*ph_support;end
        
        disp(['Phase constrained'])% range -',num2str(phase_range)])
        ph_cs=1;       %used to tell if ph constraint was used at all
        
    end
    
    if params.save_every_iterate == 1   %saves every iterate for movie making
        save_every_iterate(pnm,support,params,qq);
    end
    
    params.pn=pn;   %keep k-1
    
    %if qq <= params.charge, ALG='SF';end
    if qq == params.twin_it
        pn=next_iterate(pn,pnm,support0.*support,ALG,params.beta,params);
    else
        pn=next_iterate(pn,pnm,support,ALG,params.beta,params);
    end
    
    
    
    DM_error(qq)=calc_chi(pn0,pn);
    
    %do the averaging
    if qq >= params.iterate_avg_start %&& qq < (params.iterate_avg_start+params.iterate_avg_its)
        if params.iterate_avg_save == 1
            if calc_state(params.iterate_avg_trigger,qq) == -1
                
                avg_counter=avg_counter+1;
                disp(['Averaging estimate [',num2str(avg_counter),'/',num2str(params.iterate_avg_its),']'])
                
                sharpness(avg_counter)=sum(abs(pnm(:)).^4);  %keep 'sharpness' for selection
                
                if params.do_2D == 1,               %if it is time, keep the iterate
                   pn_g(:,:,avg_counter)=pnm; 
                else
                   pn_g(:,:,:,avg_counter)=pnm;
                end
                
                if avg_counter == numel(params.iterate_avg_trigger)
                    %[ pn_g ] = GA_breed_iterates(pn_g,sharpness,params);%
                    %old one above newer below
                    disp('Aligning all iterates (make take some time)....')
                    [ pn_g ] = align_iterates(pn_g,1,[]); 
                    pnm_avg=sum(pn_g,ndims(pn_g));
                    pn_g=[];
                end
                %pnm_avg=pnm_avg+pnm;
                %params.iterate_avg_enabled=1;  %marker to determine what to save
            end
                
        end
    end
    
            
    
end

%%

if params.GPU == 1
    try
       support0=single(support0);
    end
    data=single(data);
    support=single(support);
    pn=single(pn);
    pnm=single(pnm);
end 
params.prev_coh=[];             %can't output cell arrays to text

params.final_chi=error;
params.elapsed_time=toc;
params_out=params;
params_out.pn=0;
params_out.sw=sw;% 'Shrinkwrap was used, check shrink_trigger for details';
params_out.ph_cs=ph_cs;%'Phase constraint was used, check phase_trigger for details';

try
    if params.return_iterate == 1
        if disp_off == 0, disp('returning iterate (params.return_iterate =1)....');end
        pnm=pn;
    end
end

if ndims(pnm) < 4               %this is for simultaneous reconstructions
    params_out.coh='array';
end
    
if ndims(pnm) == 3
    coh=(fftshift(ifftn(fftshift(zero_pad_ver2(params.coh,nn(2),nn(1),nn(3))))));   % was abs
end
if ndims(pnm) == 2
    if ndims(params.coh) ==3,params.coh=extract_3D_slice(params.coh,'xy',6);end
    coh=(fftshift(ifftn(fftshift(zero_pad_ver2(params.coh,nn(2),nn(1))))));   % was abs
end

if numel(pnm_avg) > 1
    disp('Returning average....')
    pnm=pnm_avg;
end

if disp_off == 0,
    disp('///////////////////////////////////////////////////////////')
    disp(['//////',' Elapsed time for algorithm is ',num2str(params.elapsed_time),' seconds. //////'])
    disp('///////////////////////////////////////////////////////////')
end

end
%% Other functions
function params=set_param_defaults(params)

%set paramater defaults
try
    params.simultaneous_coh_det;        %this is used for doing simultaneous reconstructions
catch                                   %to determine the coherence function 
    params.simultaneous_coh_det=0;
end
try                              %backwards compatability for progressively
    params.gauss_f_trigger;      %filtering the reconstruction.  ie phase 
catch                            %low res first then high
    params.gauss_f_trigger=[-1];
end
try
    params.gauss_f_init;
    params.gauss_f_sig;
catch
    params.gauss_f_init=[0.5,.5,.5]; 
    params.gauss_f_sig=0;
end
try
   params.save_every_iterate; 
catch
   params.save_every_iterate=0; 
end
try
   params.silent;
catch
   params.silent=0;
end
try
    params.iterate_avg_start;    %when to start averaging iterates
    params.iterate_avg_its;       %number of iterates to average
    params.iterate_avg_save;       %save the average iterate
catch
    params.iterate_avg_start=200;    %when to start averaging iterates
    params.iterate_avg_its=10;       %number of iterates to average
    params.iterate_avg_save=0;       %save the average iterate
end
try 
    params.iterate_avg_inter;
catch
    params.iterate_avg_inter=2;         %interval between averages
end

try
    params.no_zero;
catch
    params.no_zero=0;
end

try
    params.GPU;
catch
    params.GPU=0;
end

try
    params.pcdi;
catch
    params.pcdi=0;
    params.coh=gauss_2D(3,3,0,0);
end
try
    params.update_its;
catch
    params.update_its=25;
end

try
    params.start_pcdi;
catch
    params.start_pcdi=30;
end
try
    params.stop_pcdi;
catch
    params.stop_pcdi=200;
end
try
    params.angle_start;
catch
    params.angle_start=30;
end

try
    params.orig;
catch
    params.orig=params.pcdi;
end
try
    params.recenter_rec;
catch
    params.recenter_rec=0;
end
try
    params.mask;
catch
    params.mask=0;
end
params.iterate_avg_trigger=0:params.iterate_avg_inter:(params.iterate_avg_inter*(params.iterate_avg_its-1));
params.iterate_avg_trigger=params.iterate_avg_trigger+params.iterate_avg_start;
try
    params.save_to_vtk;             %do the matlab to vtk conversion
catch
    params.save_to_vtk = 1;         %default is yes (not for 2d though)
end
try
    params.save_det2lab;            %do the coordinate transform
catch
    params.save_det2lab = 1;
end
try
    params.calc_resolution;
catch
    params.calc_resolution = 1;
end

try
    params.sw_min_area;
    if isempty(params.sw_min_area) == 1,params.sw_min_area=0;end
    if params.sw_min_area < 0,params.sw_min_area=0;end
catch
    params.sw_min_area=.025;
end

try
    params.sw_max_area;
    if isempty(params.sw_max_area) == 1,params.sw_max_area=1;end
    if params.sw_max_area < 0,params.sw_max_area=1;end
catch
    params.sw_max_area=.1;
end

try
    params.do_2D;
catch
    params.do_2D=0;
end


end
%%
function save_every_iterate(pnm,support,params,qq)
   
try
    
   dir_name='movie_temp/'; 
   movie_dir=[params.save_dir,dir_name]; 
   
   if qq == params.iterations,
       phasing_params_out( params,movie_dir );
       output_python_script(movie_dir,'Make_movie_reconstructing.py'); 
   end
   
   
   if isdir(movie_dir) ==0,mkdir(movie_dir);end
   
   mov_name=[num2str(qq)];
   
   while numel(mov_name) < numel(num2str(params.iterations))
       mov_name=['0',mov_name];
   end
   
   
   array=center_array(abs(pnm));
   %flip arrays for correct orientation on output
   flip=1;

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
      nf=size(array);
      nf=nf(3);

      for ff = 1:nf,        
          array(:,:,ff)=flipud(array(:,:,ff));
      end
   end
   if flop == 1
      nf=size(array);
      nf=nf(3);

      for ff = 1:nf, 
          array(:,:,ff)=flipud(rot90(array(:,:,ff),2));
      end
   end
   
   array=array/max(abs(array(:)));
   
   save([movie_dir,mov_name,'.rec'],'array'); 
   
   
   
catch
   disp('Unable to save every iterate for a movie.')
   disp('Probably cannot find params.save_dir.  Use ')
   disp('example_batch_multi_variables.m to change ')
   disp('Matlab_phasing_ver1_1.m to include the line')
   disp('params.save_dir=save_dir')
end

%if strcmp(params.ALG1,'lr') == 1 || strcmp(params.ALG2,'lr') == 1
%    
%    if exist(params,'nonGA_lres_det') ~= 1
%        
%    
%    if numel(params.nonGA_lres_det
%    
%end


end
