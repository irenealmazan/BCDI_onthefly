function [pnm support params_out chi DM_error coh] = iterative_phasing_general(pn,data,support, params)
%jclark

try
    params.GA;     
catch
    params.GA=0;     
end

if numel(params.start_guess) > 20 %doesn't have to be 20 could be anything > 1
    params.twin=0;
    disp('Setting twin=0 since a previous reconstruction will be loaded....')
end

params=set_param_defaults(params);

%% Optimise the FFTW

fftw('planner', 'patient');

%%
if params.GA ~= 1

    if params.pcdi == 1
        if ndims(data) == 2,params.coh=extract_3D_slice(params.coh,'xy');end
        [pnm support params_out chi DM_error coh ] = iterative_phasing(pn,data,support, params);
    else
        coh=zeros(size(data));          %just pass in an empty array for the coh func.  not going to be used
        params.coh=coh;             
        [pnm support params_out chi DM_error] = iterative_phasing(pn,data,support, params);
        
    end
        
elseif params.GA == 1
    
    [pn_best pn_g chi_ga params_out supports cohs]=iterative_phasing_GA(pn,data,support, params);
    
    %check what the final metric should be
    params_out.final_chi_val=chi_ga(end,:,end);
    disp('Final averaging/selection using ....')
    disp(params_out.GA_metric_return)
    switch params_out.GA_metric_return                 %decide how to choose the 'best'

        case 'sharpness'
            chi_ga=params_out.sharpness;
        case 'summed_phase'
            chi_ga=params_out.summed_phase;
        case 'area'
            chi_ga=params_out.area;
        case 'chi'
            chi_ga=chi_ga(:,:,end);

        case 'TV'
            chi_ga=params_out.TV_norm;
    end
    
    %chi_l=chi_ga(end,:,end);
    chi_l=chi_ga(end,:);
    ind=find(chi_l == min(chi_l(:)));
    
    chi=chi_ga(end,ind,:);
    chi=chi(:);
    DM_error=chi;
    
    if params.population ~= 1       %if only returning 1, don't need to worry about which to return
    
        switch params.GA_return

            case {'avg','average'}
                [summed sigma]=sum_align_GA_reconstructions(pn_g,ind);
                pnm=summed;

                if ndims(data) == 2,
                    coh=sum(cohs,3);
                    support=supports(:,:,ind);end
                if ndims(data) == 3,
                    coh=sum(cohs,4);
                    support=supports(:,:,:,ind);end

            case {'best','bst','best-cull'}

                if ndims(data) == 2,
                    pnm=pn_g(:,:,ind);
                    coh=cohs(:,:,ind);
                    support=supports(:,:,ind);end

                if ndims(data) == 3 || params.population == 1,
                    pnm=pn_g(:,:,:,ind);
                    coh=cohs(:,:,:,ind);
                    support=supports(:,:,:,ind);end

            case 'all'


            case {'avg-half','average-half','avg-half-cull'}
                avg_n=round(params.population/2);  %the number to average

                [B IX]=sort(chi_l);                %sort the chi values and get the indices

                if ndims(pn_g) == 3,
                    coh=sum(cohs(:,:,IX(1:avg_n)),3);
                    pn_g=pn_g(:,:,IX(1:avg_n));end
                
                if ndims(pn_g) == 4,                    
                    coh=sum(cohs(:,:,:,IX(1:avg_n)),4);
                    pn_g=pn_g(:,:,:,IX(1:avg_n));end

                if sum(prod(size(data))-prod(size(squeeze(pn_g)))) ~= 0
                    [summed sigma]=sum_align_GA_reconstructions(pn_g,1);
                    pnm=summed;
                else pnm=pn_g;end
                
                if ndims(data) == 2,support=supports(:,:,IX(1));end
                if ndims(data) == 3,support=supports(:,:,:,IX(1));end

            otherwise
                string=params.GA_return;
                %get the numbers from the string
                str=regexp(string,'\d','match');
                num=[];
                for ss=1:size(str,2),num=[num,char(str(ss))];end
                avg_n=str2num(num);     %thr number to average
                if avg_n > params.population
                    disp('Average number larger than population.  Averaging entire population....')
                    avg_n=params.population;
                    params_out.GA_return=['avg-',num2str(fix(avg_n))];
                else
                    disp(['Averaging ',num2str(avg_n),' iterates....'])
                end

                [B IX]=sort(chi_l);                %sort the chi values and get the indices

                if ndims(pn_g) == 3,                %2d case with several in pop
                    
                    coh=sum(cohs(:,:,IX(1:avg_n)),3);
                    support=supports(:,:,IX(1));
                    pn_g=pn_g(:,:,IX(1:avg_n));end

                if ndims(pn_g) == 4,
                    coh=sum(cohs(:,:,:,IX(1:avg_n)),4);
                    support=supports(:,:,:,IX(1));  %3d case with several in the pop
                    pn_g=pn_g(:,:,:,IX(1:avg_n));end

                if avg_n ~= 1                       %don't do if only 1 in the pop or averaging 1
                    [summed sigma]=sum_align_GA_reconstructions(pn_g,1);
                    pnm=summed;                               
                else
                    pnm=pn_g;                       %if only 1 to avg just take that 
                    
                end

        end
        params_out.coh=coh;  %keep the averaged
        
        if params.GA_save_nbest ~= 0
            file = create_save_name( params_out.files,params_out.ALG1,params_out.ALG2,params_out.iterations,params_out.pcdi,params_out.sw,params_out.seq_let,params_out.pcdi_type,params_out.GPU,params_out );
            save_dir=[params_out.save_dir,file,'/'];
            if isdir(save_dir) ==0,mkdir(save_dir);end
            save([save_dir,'pn_g.mat'],'pn_g')
        end
        
            
    else
       pnm=pn_g;                       %if only 1 to avg just take that 
       support=supports; 
       params_out.GA_return=['avg-',num2str(fix(1))];
       params_out.coh=cohs;  %will only be 1 in the cohs array
       
      
    end
    
nn=params.sz;
if ndims(data) == 2            %take the coh func to the sample plane
   coh=(fftshift(fftn(fftshift(zero_pad_ver2(params_out.coh,nn(2),nn(1))))));   %
end
if ndims(data) == 3
   coh=(fftshift(fftn(fftshift(zero_pad_ver2(params_out.coh,nn(2),nn(1),nn(3)))))); 
end

       
    
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
    params.GPU;
catch
    params.GPU=0;
end

try
    params.pcdi;
    try
        params.kernalxy;
    catch
        params.kernalxy=[0.5,0.5,0.5];
    end
    try
        params.coh;
    catch
        params.coh=gauss_3D(11,11,11,params.kernalxy(1),params.kernalxy(2),params.kernalxy(3));
    end
catch
    params.pcdi=0;
    params.coh=gauss_2D(3,3,0,0);
end

try
    params.orig;
catch
    params.orig=params.pcdi;
end

try
    params.mask;
catch
    params.mask=0;
end
params.iterate_avg_trigger=0:params.iterate_avg_inter:(params.iterate_avg_inter*(params.iterate_avg_its-1));
params.iterate_avg_trigger=params.iterate_avg_trigger+params.iterate_avg_start;


try
    params.recenter_rec;
catch
    params.recenter_rec=1;
end

%few default checks
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
    params.regularized_amp;
catch
    params.regularized_amp='none';
end

try
    params.pcdi_type;
catch
    params.pcdi_type='none';
end

try
    params.GA_save_nbest
catch
    params.GA_save_nbest=0;
end
end

