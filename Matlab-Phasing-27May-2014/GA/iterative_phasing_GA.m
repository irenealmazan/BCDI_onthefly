function [pn_best pn_g chi_ga params_out supports cohs]=iterative_phasing_GA(pn,data,support, params)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

params=set_GA_defualts(params);   %set some defaults if they are not already specified

switch lower(params.GA_return)
    
    case {'best-cull','avg-half-cull'}
        params.GA_cull_pop=1;
        
        disp('Population will be culled....')
    otherwise
        params.GA_cull_pop=0;
        disp('Population will NOT be culled....')
end
disp(' ')
disp(['Doing ',num2str(params.generations),' Generations'])
disp(['Using ',num2str(params.population),' random starts'])
disp('')
%% GA
% seed
rsd=params.GA_random_seed;

if verLessThan('matlab', '8.1.0.604')
    
    s = RandStream('mcg16807','Seed',rsd);
    RandStream.setDefaultStream(s);

else
    rng(rsd);% -- Put code to run under MATLAB 7.0.1 and later here --
end

% setup intial pn
npop=params.population;                         %n random starts
ngens=params.generations;                       %this many cycles of breeding
iterations=params.GA_iterations;                %overwrite the iterations with that for the GA
params.iterations=params.GA_iterations;         % see above comment
pn_g=zeros([params.sz,npop]);                   %create empty array for reconstructions
supports=zeros([params.sz,npop]);               %keep supports (for SW)

if ndims(data) == 3 %create array for coh fnuc and populate with initial guess. note initila guess is only used if the update is after pcdi is turned on
    cohs=zeros([params.roi(2),params.roi(1),params.roi(3),npop]); %keep the coh funcs
    for qq=1:npop,cohs(:,:,:,qq)=zero_pad_ver3(params.coh,params.roi(1),params.roi(2),params.roi(3));end
else
    cohs=zeros([params.roi(2),params.roi(1),npop]);
    for qq=1:npop,cohs(:,:,qq)=zero_pad_ver3(params.coh(:,:,6),params.roi(1),params.roi(2));end
end

chi_ga=zeros([ngens,npop,iterations]);          %keep all the chi's

params.trigger=params.GA_trigger;               %user the GA trigger, overwrite the other

L1_norm=zeros([ngens,npop]);                    %some metrics for choosing the 'best'
sharpness=zeros([ngens,npop]);                  %can easily add additional ones
summed_phase=zeros([ngens,npop]);
area=zeros([ngens,npop]);
TV_norm=zeros([ngens,npop]);
entropy=zeros([ngens,npop]);
chi_gen=zeros([ngens,npop]);

if npop == 1,params.silent=0;end                %turn screen output on if only doing 1 member of the pop

dtot=sum(data(:)).^2;
params.dtot=dtot;
params.ndims_data=ndims(data);

if ndims(data) == 3,                            %create initial starts
    
    for qq=1:npop
      if params.GA_data_rand == 0
          if qq == 1,disp('Using random support as initial guess....'),end
          pn_g(:,:,:,qq)=support.*random('uniform',.95,1,[params.sz]); 
          ptot=sum(sum(sum(abs(pn_g(:,:,:,qq)).^2)));
          pn_g(:,:,:,qq)=pn_g(:,:,:,qq)*sqrt(dtot/ptot);
      else
          
          if qq == 1,disp('Using data with random phase as initial guess....'),end
          temp=data.*exp(i*2*pi*random('uniform',0,1,[params.sz]));
          temp=ifftshift(ifftn(fftshift(temp)));
          pn_g(:,:,:,qq)=support.*temp;
          clear temp
      end
          
      supports(:,:,:,qq)=support;
   
   end
else
   
   for qq=1:npop
      if params.GA_data_rand == 0
          if qq == 1,disp('Using random support as initial guess....'),end
          pn_g(:,:,qq)=support.*random('uniform',.95,1,[params.sz]); 
          ptot=sum(sum((abs(pn_g(:,:,qq)).^2)));
          pn_g(:,:,qq)=pn_g(:,:,qq)*sqrt(dtot/ptot);
      else
          
          if qq == 1,disp('Using data with random phase as intial guess....'),end
          temp=data.*exp(i*2*pi*random('uniform',0,1,[params.sz]));
          temp=ifftshift(ifftn(fftshift(temp)));
          pn_g(:,:,qq)=support.*temp;
          clear temp
      end
      supports(:,:,qq)=support;
   end
end

%if a previous reconstruction was passed through then use this as well
if isempty(params.GA_pn_start) ~= 1 
    disp('===================================================================')
    disp('Using previous reconstruction with perturbation as intial guess....')
    params.GA_pn_start=params.GA_pn_start*sqrt(dtot/sum(abs(params.GA_pn_start(:)).^2));
    
    %align the recon with the data
    disp('Aligning previous with current data....')
    psi=fftshift(fftn(fftshift(params.GA_pn_start)));
    [h k l]=register_3d_reconstruction(abs(data),abs(psi));
    psi=((sub_pixel_shift(psi,h,k,l)));
    params.GA_pn_start=fftshift(ifftn(fftshift(psi)));
    psi=[];
    
    for qq=1:npop
        
        if ndims(pn_g) == 4,pn_g(:,:,:,qq)=pn_g(:,:,:,qq).*params.GA_pn_start;end
        if ndims(pn_g) == 3,pn_g(:,:,qq)=pn_g(:,:,qq).*params.GA_pn_start;end
   
        %if ndims(pn_g) == 4,pn_g(:,:,:,qq)=params.GA_pn_start_alpha*pn_g(:,:,:,qq)+(1-params.GA_pn_start_alpha)*params.GA_pn_start;end
        %if ndims(pn_g) == 3,pn_g(:,:,qq)=params.GA_pn_start_alpha*pn_g(:,:,qq)+(1-params.GA_pn_start_alpha)*params.GA_pn_start;end
    end
    
    disp('===================================================================')
end
%%
breed_state=1;                          %for the switching, 1 is 1, -1 is 2.
params.breed_mode=params.breed_mode1;

if params.PSO == 1
    velocity=zeros(size(pn_g));             %velocity for use with PSO
end

params.psi2=1.0;                     

pn_best=zeros(size(pn_g));              %array for keeping the best iterates 
chi_best=zeros([npop,1]);
group_best=zeros(size(support));

sharpness_best=zeros([npop,1]);
summed_phase_best=zeros([npop,1]);
area_best=zeros([npop,1]);
params.n_updates=0;                            %keep track of how many times the best was updated
pcdi_orig=params.pcdi;

if params.GA_self_seed == 1,disp('Self-seeding....'),np_order=circshift(1:npop,[0 -1]);end %make a new order so it is sequential

params_in=params;

%%
for ww=1:ngens                %do it over the generations, after each is where
                              %the 'breeding' takes place
    ref=0;
    counter=0;                          %counter to tell how many were updated
        
    if pcdi_orig == 1,disp('<< PCDI is ON >>');end
    
    if ww >= params.GA_pcdi_gen_start,      %delaying the pcdi after x generations
        params.pcdi=pcdi_orig;
        params_out.pcdi=pcdi_orig;
        params.orig=pcdi_orig;
        if pcdi_orig == 1,disp('PCDI DELAY >= CURRENT GEN'),end
    else
        if pcdi_orig == 1,disp('PCDI DELAYED....'),end
        params.pcdi=0;
        params_out.params=0;
        params.orig=0;
    end
    
    if ww > params.GA_pcdi_gen_start,params.start_pcdi=0;end      %after the first gen of pcdi we want to use pcdi mod const
    
    disp(['GENERATION - ',num2str(ww)])

    if ww >1,
        if params.twin == 1
            params.twin = 0;
            disp('Switching off twin removal for generations past the first....')
        end    
    end    %swtich the twin removal off for subsequent generations
    
    
    % this is for the data mask, to fit low res first.
    if params.GA_lres == 1
        if ww == 1,
            data0=data;
            sig0=params.sigma;
        end 
        %sigxy=(ww/params.GA_lres_genstop*(1-params.GA_lres_init)+params.GA_lres_init).^params.GA_lres_pow;
        if params.GA_lres_init == 0,params.GA_lres_init=.1;end

        alpha=calc_scale_fact(params,ww);
        params.sigma = sig0/alpha;
        
        if isnan(params.sigma),params.sigma=params.GA_sig_max;end
        
        if params.sigma <= sig0,params.sigma=sig0;end
        if params.sigma > params.GA_sig_max,params.sigma=params.GA_sig_max;end
        
        if isempty(params.GA_sig_custom) ~= 1,params.sigma=params.GA_sig_custom(ww);end
        
        if params.sigma == 0,params.threshold=0;disp('<< 0 sigma, shrinkwrap off >>');else params.threshold = params_in.threshold;end
        
        if params.threshold > 0,disp(['Adjusting SW sigma - [',num2str(params.sigma),']']),end
      
        switch params.GA_lres_type
            case 'gauss'
                gmask=make_data_mask(params,ww,data0);
            case 'gauss-many'
                gmask=make_data_mask_many(params,ww,data0); %many small gaussians
            case 'laplacian'
                gmask= make_data_mask_laplace(data);
                gmask=gmask/max(gmask(:))*(1-alpha)+alpha;
            case 'sinc'
                gmask=abs(make_data_mask_sinc(params,ww,data)).^2;
                
            case 'thresh'
                
                gmask=shrink_wrap(data0,params.GA_lres_custom(ww),0,'percent');
            
            case 'thresh-gauss'
                
                gmask=make_data_mask(params,ww,data0).*shrink_wrap(data0,params.GA_lres_custom(ww),0,'percent');
                
            case 'power'
                gmask=1;
                
        end
        
      
        data=data0.*circshift(gmask,[params.GA_lres_offset]);
        
    end
    
    % this is for the data array size
    if params.GA_lres_array_size == 1
        if ww == 1 %get the initial paramters
            init_size_fac=0.25;
            final_size=size(data0);
            init_size=init_size_fac*final_size;       %assumes arrays are even
        end
        
        sz_fac_i=( (ww-1)/(ngens-1)*(1-init_size_fac)+init_size_fac);
        
        if sz_fac_i > 1,sz_fac_i=1;end
        
        new_size=round(sz_fac_i*final_size);

        if mod(new_size(1),2) == 1,new_size(1)=new_size(1)+1;end
        if mod(new_size(2),2) == 1,new_size(2)=new_size(2)+1;end
        if ndims(data) == 3,if mod(new_size(3),2) == 1,new_size(3)=new_size(3)+1;end;end

        if ndims(data) == 2,
            pn_best=zero_pad_ver3(pn_best,new_size(2),new_size(1),size(pn_best,3));
            
            data=zero_pad_ver3(data0,new_size(2),new_size(1));end
        
        if ndims(data) == 3,
            pn_best=zero_pad_ver3(pn_best,new_size(2),new_size(1),new_size(3),size(pn_best,4));
            data=zero_pad_ver3(data0,new_size(2),new_size(1),new_size(3));end
        
        params.sz=size(data);
        
        disp(['New size [x,y] -> [',num2str(size(data,2)),',',num2str(size(data,1)),']'])
        
    end
    
    for qq=1:npop                           %iterate on each of the population

        if params.GA_self_seed == 1,qqq=np_order(qq);else qqq=qq;end
        
        if ndims(data) == 3,
            support=supports(:,:,:,qq);     %use the last support
            pn=pn_g(:,:,:,qq);              %get the next iterate
            params.coh=cohs(:,:,:,qq);end   %get the coh
        if ndims(data) == 2,
            support=supports(:,:,qq);
            pn=pn_g(:,:,qq);
            params.coh=cohs(:,:,qq);end

%       %resize pn and support
        if params.GA_lres_array_size == 1
            support=shrink_wrap(ResizeFFt(zeros(size(data)),support),.2,0);            
            pn=ResizeFFt(zeros(size(data)),pn);
        end
        
        %do the iterations
        [pnm support params_out chi DM_error coh] = iterative_phasing(pn,data,support, params);

        %resize pn and support
        if params.GA_lres_array_size == 1,
            support=shrink_wrap(ResizeFFt(zeros(size(data0)),support),.2,0);
            pnm=ResizeFFt(zeros(size(data0)),pnm);
        end
        
        if params.pre_align == 1
            if qq == 1,ref=pnm;end
            pnm=align_and_zero(ref,pnm); 
        end
        
        if ndims(data) == 3,
            supports(:,:,:,qqq)=support;              %save new support
            cohs(:,:,:,qqq)=real(zero_pad_ver3(fftshift(fftn(fftshift(coh))),params.roi(1),params.roi(2),params.roi(3)));           %save updated coh func
            if ww ~= ngens
                if params.S_bfr_breed == 1
                    pn_g(:,:,:,qqq)=pnm.*support;
                else
                    pn_g(:,:,:,qqq)=pnm;
                end
            else          %save new estimate
                pn_g(:,:,:,qqq)=pnm;
            end
        end
        if ndims(data) == 2,
            supports(:,:,qqq)=support;
            cohs(:,:,qqq)=zero_pad_ver3(fftshift(ifftn(fftshift(coh))),params.roi(1),params.roi(2));          %save updated coh fu
            if ww ~= ngens
                if params.S_bfr_breed == 1
                    pn_g(:,:,qqq)=pnm.*support;
                else
                    pn_g(:,:,qqq)=pnm;
                end
            else          %save new estimate
                pn_g(:,:,qqq)=pnm;
            end
        end
        chi_ga(ww,qqq,:)=chi;                         % save chi
        
        [ out ] = calculate_things(pnm);
        
        L1_norm(ww,qqq)=out.L1_norm;
        sharpness(ww,qqq)=out.sharpness;
        summed_phase(ww,qqq)=-1*out.summed_phase; %-ve so that the largest is the minimum of the metric
        area(ww,qqq)=-1*out.area;
        chi_gen(ww,qqq)=chi(end);
        
        if strcmp(params.GA_metric,'TV')
            TV_temp=calc_TV(abs(pnm));
            TV_norm(ww,qqq)=sum(TV_temp(:));
            TV_temp=[];
        end
            
        if ww ~= 1,
            % check if the current is better, if then update best
            %if sharpness(ww,qq) < sharpness_best(qq),
            if chi_ga(ww,qqq,end) < chi_best(qqq),
                
                disp(['------UPDATING BEST [',num2str(qqq),'] --- ',num2str(chi_best(qqq)),' --> ',num2str(chi_ga(ww,qqq,end))])
                
                % update chi best
                summed_phase_best(qqq)=summed_phase(ww,qqq);
                area_best(qqq)=area(ww,qqq);
                sharpness_best(qqq)=sharpness(ww,qqq);
                chi_best(qqq)=chi_ga(ww,qqq,end);
                
                if params.GA_lres_array_size ~= 1 %don't update best if using smaller aray size, need to fix this
                    if ndims(data) == 2,pn_best(:,:,qq)=pn_g(:,:,qq);end
                    if ndims(data) == 3,pn_best(:,:,:,qq)=pn_g(:,:,:,qq);end
                end
                
                counter=counter+1;        %if any are better inc counter
                params.n_updates=params.n_updates+1;  %keep track of total
            end
        end
        
    end                 %end of the population iterations

    % after this is where the intermixing takes place
    if ww ~= 1
        if counter == 0,                %will be 0 if none were updated
            disp('NONE UPDATED')
            
            if params.GA == 1
               if params.GA_switch == 1 
                   breed_state=breed_state*(-1);     %switch if nothing is happening
                   disp('SWITCHING BREED MODES')
                   if breed_state == 1, 
                    
                       params.breed_mode=params.breed_mode1;
                       disp(['USING BREED MODE 1 - [',params.breed_mode1,']'])
                   
                   end
                   if breed_state == -1, 
                      
                       disp(['USING BREED MODE 2 - [',params.breed_mode2,']'])
                       params.breed_mode=params.breed_mode2;
                   
                   end
               end
            end
            
            if params.PSO == 1
                if params.PSO_inc_psi == 1
                    parama.psi2=params.psi2+1.0;
                    disp(['INCREASING PSI2 - [',num2str(params.psi2),']'])
                end
            end
            
        end
    end
    
    % set the first iterate to the best since we don't have one yet
    if ww == 1,
        sharpness_best=sharpness(ww,:);
        chi_best(:)=chi_ga(ww,:,end);
        pn_best=pn_g;
        disp('------SETTING FIRST ITERATE TO BEST------')
    end
    
    disp('       Best      Current')
    for vv=1:npop,disp([squeeze(chi_best(vv)),chi_ga(ww,vv,end)]),end %'[',num2str(vv),'] ',
        
    
    % find the group best
    %ind=find(chi_best == min(chi_best(:)));
    chi_temp=get_chi_choose(sharpness,summed_phase,area,chi_ga,TV_norm,ww,params);
    ind=find(chi_temp == min(chi_temp(:))); 
    
    % chi best and pn_best will have the group best information and iterate
    if ndims(data) == 2,group_best=pn_best(:,:,ind);end
    if ndims(data) == 3,group_best=pn_best(:,:,:,ind);end
    disp(['------GROUP BEST IS [',num2str(ind),']------'])
    
    % PSO
    if params.PSO == 1
        if ww ~= ngens
            disp('------PSO UPDATING ITERATES------')
            if params.use_sharp == 1,chi_choose=sharpness(ww,:,end);else
                chi_choose=chi_ga(ww,:,end);end %chi_best
            [pn_g velocity]=PSO_iterates(pn_g,pn_best,velocity,group_best,chi_choose,params);
        end
    end
   
    % GA
    if params.GA == 1
        
        %select the metric
        disp(params.GA_metric)
        chi_choose = get_chi_choose(sharpness,summed_phase,area,chi_ga,TV_norm,ww,params);
                
        %do the breeding
        if strcmp(params.breed_mode1,'none') == 0
            if ww ~= ngens
                disp('------GA UPDATING ITERATES------')
                
                if params.GA_cull_pop == 1
                    nremove=params.GA_cull_remove;      %the number to remove from the population
                    
                    if (npop-nremove) > params.npop_stop  %don't remove some if the population is too small
                                                
                        [sorted_chi pop_rank]=sort(chi_choose);    %rank them according to their fitness
                        remove_ind=pop_rank(end-nremove+1:end);
                        
                        disp(['Eliminating ',num2str(nremove),' iterates from the population....'])
                        
                        chi_choose(remove_ind)=[];

                        if ndims(data) == 2,     
                            pn_g(:,:,remove_ind)=[];
                            pn_best(:,:,remove_ind)=[];                                
                        end

                        if ndims(data) == 3,                               
                            pn_g(:,:,:,remove_ind)=[];
                            pn_best(:,:,:,remove_ind)=[];                               
                        end

                        L1_norm(:,remove_ind)=[];
                        TV_norm(:,remove_ind)=[];
                        sharpness(:,remove_ind)=[];
                        summed_phase(:,remove_ind)=[];
                        area(:,remove_ind)=[];
                        chi_ga(:,remove_ind,:)=[];
                        chi_best(remove_ind)=[];
                        chi_gen(:,remove_ind)=[];
                                          
                        npop=npop-nremove;
                        
                        disp(['New population size ',num2str(npop)])
                        
                    end
                end                
                
                if npop ~= 1 %only do the 'breeding' if there is more than 1            
                    [ pn_g ] = GA_breed_iterates_ver2(pn_g,chi_choose,params,group_best,data);%chi_ga(ww,:,end)
                end
                
            end
        else
            chi_choose=chi_ga(ww,:,end);
        end
        % save the best iterate ina temp folder
        if params.GA_save_gen_temp == 1
            output_iterate(pn_g,supports,chi_choose,params,data,chi_gen,sharpness,summed_phase,L1_norm,area,ww)
        end
        
        %this step is important for reconstructions using a shrinkwrap support and/or one that has twin removal                                   
        if params.GA_threshold ~= 0 || params.threshold ~= 0 %need to update the supports after shifting etc
            
            if params.threshold ~= 0       %use the same parameters as for the normal sw
                supports=generate_new_supports(pn_g,params.threshold,params.sigma);
            else
                supports=generate_new_supports(pn_g,params.GA_threshold,params.GA_sig); %otherwise just do it after breeding
            end
        end
        
    end
    
    
end

chi=chi_ga(end,:,end);
ind=find(chi == min(chi(:)));

pn=pn_g(:,:,ind);

params_out.sharpness=sharpness/max(sharpness(:));
params_out.summed_phase=summed_phase;
params_out.area=area;
params_out.chi=chi_ga(:,:,end);
params_out.TV_norm=TV_norm/max(TV_norm(:));



end

function params=set_GA_defualts(params)

try
    params.S_bfr_breed;
catch
    params.S_bfr_breed=0;
end

try
    params.PSO;
catch
    params.PSO = 0;
end
try
    params.pre_align;
catch
    params.pre_align=0;
end

try
    params.rand_x_only;
catch
    params.rand_x_only=0;
end
try
    params.rand_y_only;
catch
    params.rand_y_only=0;
end
%%
try
    params.GA_lres;
    params.GA_lres_init;
    params.GA_lres_genstop;
catch
    params.GA_lres=0;             %start with low res data
    params.GA_lres_init=.1;   %sig of gauss mask of array size, 0-1
    params.GA_lres_genstop=10;    %at which gerneation the data should be fill resn
end
try
    params.GA_lres_pow;         
catch
    params.GA_lres_pow=1;         %the power of the change of sig
end

try
    params.GA_metric;
catch
    params.GA_metric='sharpness';
end

try
    params.GA_metric_return;
catch
    params.GA_metric_return='chi';
end

try
    params.use_sharp;  %so it is still compatable with older versions
    params.GA_metric='sharpness';
end

try
    params.GA_threshold;                %to do a support after each generation, based on the best
catch
    params.GA_threshold=0.0;                
end
 
try
    params.GA_sig;                %to do a support after each generation, based on the best
catch
    params.GA_sig=0.5;                
end

try
    params.GA_return;
catch
    params.GA_return='avg';
end

try
    params.GA_data_rand;
catch
    params.GA_data_rand=0;
end

try
    params.npop_stop;
catch
    params.npop_stop=6;
end

try
    params.GA_lres_array_size;
catch    
    params.GA_lres_array_size=0;
end

try
    params.GA_lres_offset;
catch 
    params.GA_lres_offset=0;
end

try
    params.GA_cull_remove;
catch
    params.GA_cull_remove=2;    
end

try
    params.GA_sig_max;
catch
    params.GA_sig_max=4;
end

try
    params.GA_pcdi_gen_start;
catch
    params.GA_pcdi_gen_start=1;
end
try
    params.GA_lres_type;
catch
    params.GA_lres_type='gauss';
end

try
    params.GA_random_seed;
catch
    params.GA_random_seed=1;
end
try
    params.GA_keep_lres_phase;
catch
    params.GA_keep_lres_phase=0;    
end

try
    params.GA_sig_custom;  %custom sig values, leave empty for automatic mode
    if numel(params.GA_sig_custom) ~= params.generations,params.GA_sig_custom=[params.GA_sig_custom,(zeros(1,(params.generations))+params.GA_sig_custom(end))];end
catch
    params.GA_sig_custom=[];  %
end

try
    params.GA_lres_custom;
    if numel(params.GA_lres_custom) ~= params.generations,params.GA_lres_custom=[params.GA_lres_custom,(zeros(1,(params.generations))+params.GA_lres_custom(end))];end
catch
    params.GA_lres_custom=[];
end

try
    params.GA_self_seed;
catch
    params.GA_self_seed=0;
end

try
    params.GA_pn_start;
catch
    params.GA_pn_start=[];
end

try
    params.GA_pn_start_alpha;
catch
    params.GA_pn_start_alpha=0.1;
end

try
    params.GA_save_gen_temp;
catch
    params.GA_save_gen_temp=1;
end

try
    params.GA_save_nbest;
catch
    params.GA_save_nbest=0;
end

params.silent=1;

end

function gmask = make_data_mask(params,ww,data)


%sigxy=( (ww-1)/(params.GA_lres_genstop-1)*(1-params.GA_lres_init)+params.GA_lres_init).^params.GA_lres_pow;                
%if sigxy < params.GA_lres_init,sigxy=params.GA_lres_init;end
sigxy = calc_scale_fact(params,ww);

if isempty(params.GA_lres_custom) ~= 1,sigxy=params.GA_lres_custom(ww);disp('Using custom lres values....');end

if sigxy < 1,
    disp(['USING GAUSS MASK - ',num2str(sigxy)])
    if ndims(data) == 2,
        gx=size(data,2);
        gy=size(data,1);
        gmask=gauss_2D(gx,gy,gx*sigxy,gy*sigxy);
        %gmask=sinc_2D(gx,gy,sigxy*2,sigxy*2);
        %gmask=super_gauss_2D(gx,gy,gx*sigxy,gy*sigxy,2);
    end

    if ndims(data) == 3,
        gx=size(data,2);
        gy=size(data,1);
        gz=size(data,3);
        gmask=gauss_3D(gx,gy,gz,gx*sigxy,gy*sigxy,gz*sigxy);
    end
else
    gmask=1; 
end


end

function gmask = make_data_mask_many(params,ww,data)

params=set_mask_defaults(params);  %set defaults if non set

sigxy=params.GA_lres_cust_sigxy;

rad=ww/params.GA_lres_genstop*max(size(data));

ngpr=1.0/sigxy; %number per row
ngauss=(ngpr)^(ndims(data));  %number to use

%make basis gaussian
if ndims(data) == 2,
    gx=size(data,2);
    gy=size(data,1);
    gz=size(data,3);
    gbasis=gauss_2D(gx,gy,gx*sigxy,gy*sigxy);
end
if ndims(data) == 3,
    gx=size(data,2);
    gy=size(data,1);
    gz=size(data,3);
    gbasis=gauss_3D(gx,gy,gz,gx*sigxy,gy*sigxy,gz*sigxy);
end

lattice=create_nd_lattice(ngpr,data).*generate_circle_nd(size(data),rad,ndims(data));
gmask=zeros(size(lattice));
ind=find(lattice > 0);
[x y z]=ind2sub(size(lattice),ind);

for et=1:numel(x),gmask=gmask+circshift(gbasis,[round(gy/2)-y(et),round(gx/2)-x(et),round(gz/2)-z(et)]);end
gmask=gmask/max(gmask(:));

end

function params=set_mask_defaults(params)

try
    params.GA_lres_cust_sigxy;
catch
    params.GA_lres_cust_sigxy=.05;  %fraction of array size
end



end

function gmask = make_data_mask_sinc(params,ww,data)

sigxy = calc_scale_fact(params,ww);

if sigxy < 1,
    disp(['USING GAUSS MASK - ',num2str(sigxy)])
    if ndims(data) == 2,
        gx=size(data,2);
        gy=size(data,1);
        
        [x,y]=meshgrid( -(gx-1)/2:(gx-1)/2,-(gy-1)/2:(gy-1)/2);
        x=x/max(x(:))*.5;
        y=y/max(y(:))*.5;
        gmask=sinc( sqrt(abs(x).^2+abs(y).^2)*(1-sigxy)*.5*max([gx,gy]) );
        
    end

    if ndims(data) == 3,
        gx=size(data,2);
        gy=size(data,1);
        gz=size(data,3);
        gmask=gauss_3D(gx,gy,gz,gx*sigxy,gy*sigxy,gz*sigxy);
    end
else
    gmask=1; 
end


end

function gmask = make_data_mask_lin(params,ww,data)

sigxy = calc_scale_fact(params,ww);

if sigxy < 1,
    disp(['USING GAUSS MASK - ',num2str(sigxy)])
    if ndims(data) == 2,
        gx=size(data,2);
        gy=size(data,1);
        gmask=gauss_2D_v2(gx,gy,gx,gy,sigxy);
    end

    if ndims(data) == 3,
        gx=size(data,2);
        gy=size(data,1);
        gz=size(data,3);
        gmask=gauss_3D(gx,gy,gz,gx,gy,gz,sigxy);
    end
else
    gmask=1; 
end


end

function sigxy = calc_scale_fact(params,ww)

sigxy=( (ww-1)/(params.GA_lres_genstop-1)*(1-params.GA_lres_init)+params.GA_lres_init).^params.GA_lres_pow;
                
if sigxy < params.GA_lres_init,sigxy=params.GA_lres_init;end

end

function kern = make_data_mask_laplace(data)

nx=size(data,2);
ny=size(data,1);
if ndims(data) == 3,nz=size(data,3);end

h = zeros(3,3,3);
h(:,:,1) = [0 3 0;3 10 3;0 3 0];
h(:,:,3) = h(:,:,1);
h(:,:,2) = [3 10 3;10 -96 10;3 10 3];

if ndims(data) == 3,ff=zero_pad_ver3(h,nx,ny,nz);end
if ndims(data) == 2,ff=zero_pad_ver3(h(:,:,2),nx,ny);end

kern=abs(fftshift(fftn(fftshift(ff/sum(abs(ff(:)))))));

end

function output_iterate(pn_g,supports,chi_choose,params,data,chi_gen,sharpness,summed_phase,L1_norm,area,ww)


if isdir([params.save_dir,'temp/']) == 0,mkdir([params.save_dir,'temp/']);end

ind=find(chi_choose == min(chi_choose));

if ndims(data) == 3,
    best=pn_g(:,:,:,ind);
    ss=supports(:,:,:,ind);  
end

if ndims(data) == 2,
    best=pn_g(:,:,ind);
    ss=supports(:,:,ind);  
end

[best yxz]=center_array(best);
ss=circshift(ss,yxz);

if ndims(data) == 3
    save_display_rec3D(best,ss,[params.save_dir,'temp/',num2str(ww)] );
end
if ndims(data) == 2
    save_display_rec2D(best,ss,[params.save_dir,'temp/',num2str(ww)] );
end

save([params.save_dir,'temp/chi_gen.mat'],'chi_gen')
save([params.save_dir,'temp/sharpness.mat'],'sharpness')
save([params.save_dir,'temp/summed_phase.mat'],'summed_phase')
save([params.save_dir,'temp/L1_norm.mat'],'L1_norm')
save([params.save_dir,'temp/area.mat'],'area')

end

function chi_choose = get_chi_choose(sharpness,summed_phase,area,chi_ga,TV_norm,ww,params)

switch params.GA_metric                 %decide how to choose the 'best'

    case 'sharpness'
        chi_choose=sharpness(ww,:);
    case 'summed_phase'
        chi_choose=summed_phase(ww,:);
    case 'area'
        chi_choose=area(ww,:);
    case 'chi'
        chi_choose=chi_ga(ww,:,end);

    case 'TV'
        chi_choose=TV_norm(ww,:);

end

end

function fxy = gauss_2D_v2(nx,ny,sigx,sigy,alpha)
%
try
    alpha;
catch
    alpha=1;
end

fxy=zeros(nx,ny);

[x , y]=meshgrid( -(nx-1)/2:(nx-1)/2,-(ny-1)/2:(ny-1)/2);

if sigx < 0
    if mod(nx,2) == 1
        if sigx == 0, gx=(x == 0);else gx=exp(alpha*0.5.*x.^2./sigx^2);end
    else gx=exp(alpha*0.5.*x.^2./sigx^2);end
else
    if mod(nx,2) == 1
        if sigx == 0, gx=(x == 0);else gx=exp(-alpha*0.5.*x.^2./sigx^2);end
    else gx=exp(-alpha*0.5.*x.^2./sigx^2);end
end

if sigy < 0
    if mod(ny,2) == 1
        if sigy == 0, gy=(y == 0);else gy=exp(alpha*0.5.*y.^2./sigy^2);end
    else gy=exp(alpha*0.5.*y.^2./sigy^2);end
else
    if mod(ny,2) == 1
        if sigy == 0, gy=(y == 0);else gy=exp(-alpha*0.5.*y.^2./sigy^2);end
    else gy=exp(-alpha*0.5.*y.^2./sigy^2);end    
end

fxy=gx.*gy;

fxy=fxy/sum((fxy(:)));

end

