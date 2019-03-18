function [ coh, params ] = determine_coh(Mk,Mm,params )
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

coh_dim=ndims(Mk);      %the dimensions of the final coh func

gn=11 ;      %size of gaussian kernal

params = set_determine_coh_defualts(params);

if params.pcdi_seperable == 1, disp('Assuming coh func is seperable....'),end
%xyz=params.sz/2-params.roi;

%% do only required dimensions
if numel(params.do) == 2, params.do=[params.do,1];end    %set z = not do if its 2d

if sum(params.do) == 0          % do them all if specified

    if ndims(Mm) == 2   
        if sum(params.roi(:) <= 0)
            disp('Using full ROI....')
            Mm_r=Mm;
            Mk_r=Mk;
        else
            Mm_r=zero_pad(Mm,params.roi(1),params.roi(2)); %replaced extract_max with extract_center
            Mk_r=zero_pad(Mk,params.roi(1),params.roi(2));
        end
    end
    if ndims(Mm) == 3    
        if sum(params.roi(:) <= 0)
            disp('Using full ROI....')
            Mm_r=Mm;
            Mk_r=Mk;
        else
            Mm_r=zero_pad(Mm,params.roi(1),params.roi(2),params.roi(3));

            try
                Mk_r=zero_pad(Mk,params.roi(1),params.roi(2),params.roi(3));
            catch
                disp('- Warning - Warning - Warning - Warning - Warning - Warning - Warning -')
                disp('Iterates detector wavefield has a maximum that is not centered....')
                disp('Assuming that the centre of the wavefield is the array center....')
                disp('- Warning - Warning - Warning - Warning - Warning - Warning - Warning -')
                size_Mk=size(Mk);
                lind=round(size_Mk/2)-params.roi/2;
                hind=round(size_Mk/2)+params.roi/2;
                Mk_r=Mk(lind(1):hind(1),lind(2):hind(2),lind(3):hind(3));
                %Mk_r=extract_max(temp,params.roi(1),params.roi(2),params.roi(3));
                %temp=0;
                lind=0;
                hind=0;
            end
        end
        
    end
        
end


if sum(params.do) ~= 0          % if any are specified to be not done do 2D
    
    sz=size(Mm);
    
    if max(size(sz(:))) == 2, sz=[sz,1.1];end   %check for 2d, make sz >1 so that round(sz/2)=1
    
    if sum(params.do) == 1
        nd=find(params.do); % find the non-zero element, i.e the not one
        if nd == 1, 
            Mm=squeeze(Mm(:,round(sz(2)/2),:));%zy
            Mk=squeeze(Mk(:,round(sz(2)/2),:));%zy
            orien0='zy';
            Mm_r=zero_pad(Mm,params.roi(3),params.roi(2));
            Mk_r=zero_pad(Mk,params.roi(3),params.roi(2));
        end  
        if nd == 2,
            %Mm=squeeze(Mm(round(sz(1)/2),:,:));%zx
            %Mk=squeeze(Mk(round(sz(1)/2),:,:));%zx
            %orien0='zx';
            
            Mm=rot90(squeeze(Mm(round(sz(1)/2),:,:)));%xz
            Mk=rot90(squeeze(Mk(round(sz(1)/2),:,:)));%xz
            orien0='xz';
            
            Mm_r=zero_pad(Mm,params.roi(3),params.roi(1));
            Mk_r=zero_pad(Mk,params.roi(3),params.roi(1));
        end  
        if nd == 3, 
            Mm=squeeze(Mm(:,:,round(sz(3)/2)));%xy
            Mk=squeeze(Mk(:,:,round(sz(3)/2)));%xy
            orien0='xy';
            Mm_r=zero_pad(Mm,params.roi(1),params.roi(2));
            Mk_r=zero_pad(Mk,params.roi(1),params.roi(2));
        end  
        do=[0,0,0]; %since it is going to be 2d anyway, set them all to 0
    end
    
    if sum(params.do) >= 2
        nd=find(params.do == 0);    %find the one to do
        if nd == 1, 
            
            if params.pcdi_seperable == 1
                Mm=sum(Mm,3);
                Mk=sum(Mk,3);
            else
                Mm=squeeze(Mm(:,:,round(sz(3)/2)));
                Mk=squeeze(Mk(:,:,round(sz(3)/2)));
            end
            
            orien0='xy';
            do=[0,1,1];
            Mm_r=zero_pad(Mm,params.roi(1),params.roi(2));
            Mk_r=zero_pad(Mk,params.roi(1),params.roi(2));
        end  %xy , do x 
        if nd == 2, 
            
            if params.pcdi_seperable == 1
                Mm=sum(Mm,3);
                Mk=sum(Mk,3);
            else
                Mm=squeeze(Mm(:,:,round(sz(3)/2)));
                Mk=squeeze(Mk(:,:,round(sz(3)/2)));
            end
                
            orien0='xy';
            do=[1,0,1];
            Mm_r=zero_pad(Mm,params.roi(1),params.roi(2));
            Mk_r=zero_pad(Mk,params.roi(1),params.roi(2));
        end  %xy , do y
        if nd == 3, 
            %Mm=squeeze(Mm(round(sz(1)/2),:,:)); %zx
            %Mk=squeeze(Mk(round(sz(1)/2),:,:)); %zx
            %orien0='zx';
            %do=[0,1,1];
            if params.pcdi_seperable == 1
                Mm=rot90(squeeze(sum(Mm,1)));%xz
                Mk=rot90(squeeze(sum(Mk,1)));%xz
            else
                Mm=rot90(squeeze(Mm(round(sz(1)/2),:,:)));%xz
                Mk=rot90(squeeze(Mk(round(sz(1)/2),:,:)));%xz
            end
                
            orien0='xz';
            do=[1,0,1];
            
            Mm_r=zero_pad(Mm,params.roi(3),params.roi(1));
            Mk_r=zero_pad(Mk,params.roi(3),params.roi(1));
        end   %zx , do z
    end
    
    lb=determine_sigmas_in(orien0,params.min_sigmas);
    ub=determine_sigmas_in(orien0,params.max_sigmas);

else                %doing all 3 dims
    do=params.do;
    
    if ndims(Mm) == 3
        %Mm_r=zero_pad(Mm,params.roi(1),params.roi(2),params.roi(3));
        %Mk_r=zero_pad(Mk,params.roi(1),params.roi(2),params.roi(3));
        lb=params.min_sigmas;
        ub=params.max_sigmas;
    end
    if ndims(Mm) == 2
        %Mm_r=zero_pad(Mm,params.roi(1),params.roi(2));
        %Mk_r=zero_pad(Mk,params.roi(1),params.roi(2));
        orien0='xy';
        lb=determine_sigmas_in(orien0,params.min_sigmas);
        ub=determine_sigmas_in(orien0,params.max_sigmas);
    end
end



if params.pcdi_norm == 1,
    disp('Normalising...')
    Mk_r=sqrt(Mk_r.^2/sum(sum(sum(Mk_r.^2)))*sum(sum(sum(Mm_r.^2))));
end
%###################LUCY DECONV
%% Lucy deconv
lucy_its=params.lucy_its;

try
    kern_size=params.conv_kern_size;
catch
    kern_size=17;
end

try
    params.prev_coh;
catch
    params.prev_coh={Mm_r.^2};
end

if params.symmetrize_data == 1
   disp('Symmetrizing data for symmetrical coherence function....')
   Mm_r=center_array(Mm_r)+center_array(flipdim(flipdim(flipdim(Mm_r,1),2),3));
   Mk_r=center_array(Mk_r)+center_array(flipdim(flipdim(flipdim(Mk_r,1),2),3));
end

if strcmp(params.pcdi_type,'lucy')
    
    disp(['Updating coherence function using iterative lucy-richardson',' [',num2str(lucy_its),' iterations]...'])
    
    if params.use_previous_coh == 0
        coh=deconvlucy(Mm_r.^2,Mk_r.^2,lucy_its);
    else
        disp('Using previous RL output as input....')
        params.prev_coh=deconvlucy(params.prev_coh,Mk_r.^2,lucy_its);
        coh=cell2mat(params.prev_coh(2));
    end
        
    if ndims(coh) == 2
        coh=center_array(coh);
        
        %coh=extract_max(coh,kern_size,kern_size);
        if kern_size >=0,coh=extract_max(coh,kern_size,kern_size);else
            disp('Kern size = ROI size....')
        end
        %coh2d=symmetrize_kernal(coh);
        coh2d=coh;%
        
        if coh_dim == 3,coh=make_slice_3D(coh2d,zeros([kern_size,kern_size,kern_size]),orien0,(kern_size+1)/2);end
    
    end
    
    %if params.use_fftconv == 0,
    if ndims(coh) == 3,
        if kern_size >=0,coh=extract_max(coh,kern_size,kern_size,kern_size);else
            disp('Kern size = ROI size....'),
            %coh=center_array(coh);      %just to be sure its centered!
        end
    end
    %end
    
    if params.symmetrize_kernal == 1,
        disp('Symmetrizing coherence function....')
        coh=symmetrize_kernal(coh);
    end
    
    coh=abs(coh)/sum(sum(sum(abs(coh))));
    params.coh=coh;
    
end

% Lucy deconv with additional minimization step, use pattern search
if strcmp(params.pcdi_type,'lucy_ps')
    disp('Updating coherence function using iterative lucy-richardson with additional refinement from pattern searching...')
    
    %do the lucy deconv intialliy
    coh=deconvlucy(Mm_r.^2,Mk_r.^2,lucy_its);
    
    if ndims(coh) == 2
        coh=extract_max(coh,7,7);
        coh=abs(coh)/sum(sum(abs(coh)));
        coh2d=coh;
    else
        coh=extract_max(coh,7,7,7);
        coh=abs(coh)/sum(sum(sum(abs(coh))));
    end
    %now update the values using a direct minimisation scheme
    options=psoptimset('Display','off','TolFun',1e-7,'TolMesh',1e-6,'MeshContraction',.5,'PollingOrder','Random','MeshExpansion',3);
    x0=coh;
    lb=zeros(size(x0));
    ub=ones(size(x0));
    f=@(x)objfun_kern2D(x,Mm_r,Mk_r);               %can actually do 3d as well
    x=patternsearch(f,x0,[],[],[],[],lb,ub,[],options);
    
    if ndims(coh) == 2
        coh2d=padarray(x,[2,2]);
        coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    else
        coh=padarray(x,[2,2,2]);
    end
    params.coh=coh;
    
end
if strcmp(params.pcdi_type,'lucy_lsqs')
    disp('Updating coherence function using iterative lucy-richardson with additional refinement from least squares...')
    
    %do the lucy deconv intialliy
    coh=deconvlucy(Mm_r.^2,Mk_r.^2,500);
    
    if ndims(coh) == 2
        coh=extract_max(coh,7,7);
        coh=abs(coh)/sum(sum(abs(coh)));
        coh2d=coh;
    else
        coh=extract_max(coh,7,7,7);
        coh=abs(coh)/sum(sum(sum(abs(coh))));
    end
    %x0=extract_3D_slice(params.coh,orien0,6 );
    x0=coh;
    options = optimset('TolFun',1e-10,'Algorithm','trust-region-reflective','TolX',1e-8);
    lb=zeros(size(x0));
    ub=ones(size(x0));
    
    f=@(x)least_squares(x,Mm_r,Mk_r);
    [x resnorm residual]=lsqnonlin(f,x0,lb,ub,options);
    
    if ndims(coh) == 2
        coh2d=x;
        coh2d=symmetrize_kernal(coh2d);
        coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    else
       coh=x;
    end
    coh=abs(coh)/sum(sum(sum(abs(coh))));    
    params.coh=coh;
    
end
%###################Weiner DECONV - relies on using previous (k-1)
if strcmp(params.pcdi_type,'weiner') || strcmp(params.pcdi_type,'trunc-weiner')
    
    things=32;  %fourier 'modes' to keep    
        
    if strcmp(params.pcdi_type,'weiner'),
        disp(['Updating coherence function using Weiner filter.....'])
    end
    if strcmp(params.pcdi_type,'trunc-weiner'),
        disp(['Updating coherence function using truncated Weiner filter.....'])
    end
    
    coh=deconvwnr(Mm_r.^2,Mk_r.^2);
    coh(coh < 0)=0;
        
    if strcmp(params.pcdi_type,'trunc-weiner'),
        sorted_vals=sort(coh(:),'descend');
        coh( coh < sorted_vals(things))=0;
    end
    
    if ndims(coh) == 2
        coh=center_array(coh);
        coh=extract_max(coh,kern_size,kern_size);
        %coh2d=symmetrize_kernal(coh);
        coh2d=coh;%
        
        if coh_dim == 3,coh=make_slice_3D(coh2d,zeros([kern_size,kern_size,kern_size]),orien0,(kern_size+1)/2);end
    
    end
    
    %if params.use_fftconv == 0,
    if ndims(coh) == 3,
        if kern_size >=0,coh=extract_max(coh,kern_size,kern_size,kern_size);else
            disp('Kern size = ROI size....'),end
    end
    %end
    
    if params.symmetrize_kernal == 1,
        disp('Symmetrizing coherence function....')
        coh=symmetrize_kernal(coh);
    end
    
    coh=(coh)/sum(sum(sum(abs(coh))));
    params.coh=coh;
    
end
%#################### Weiner and RL deconv, seems to work for 2k-1
if strcmp(params.pcdi_type,'lucy-weiner')
    
    disp(['Updating coherence function using Weiner filter & iterative RL [',num2str(lucy_its),']'])
    
    
    cohW=deconvwnr(Mm_r.^2,Mk_r.^2);
    cohW(cohW < 0)=0;
    
    if params.use_previous_coh == 0
        coh=deconvlucy(Mm_r.^2,Mk_r.^2,lucy_its);
    else
        disp('Using previous RL output as input....')
        params.prev_coh=deconvlucy(params.prev_coh,Mk_r.^2,lucy_its);
        coh=cell2mat(params.prev_coh(2));
    end
    
    coh=abs(coh.*cohW).^.5;
        
    if ndims(coh) == 2
        coh=center_array(coh);
        coh=extract_max(coh,kern_size,kern_size);
        %coh2d=symmetrize_kernal(coh);
        coh2d=coh;%
        
        if coh_dim == 3,coh=make_slice_3D(coh2d,zeros([kern_size,kern_size,kern_size]),orien0,(kern_size+1)/2);end
    
    end
    
    %if params.use_fftconv == 0,
    if ndims(coh) == 3,
        if kern_size >=0,coh=extract_max(coh,kern_size,kern_size,kern_size);else
            disp('Kern size = ROI size....'),end
    end
    %end
    
    if params.symmetrize_kernal == 1,
        disp('Symmetrizing coherence function....')
        coh=symmetrize_kernal(coh);
    end
    
    coh=abs(coh)/sum(sum(sum(abs(coh))));
    params.coh=coh;
    
end
%%##
if strcmp(params.pcdi_type,'HyBR')
    
    disp(['Updating coherence function using HyBR...'])
    
    PSF=Mk_r.^2;
    IPC=Mm_r.^2;
    %   Use RestoreTools to construce matrix and preconditioner objects:
    A = psfMatrix(PSF,'zero');
    P = psfPrec(A, IPC);
    %
    %   Use HyBR with default settings and preconditioning to solve
    [coh, output] = HyBR(A, IPC, P);
    
        
    if ndims(coh) == 2
        coh=center_array(coh);
        coh=extract_max(coh,kern_size,kern_size);
        %coh2d=symmetrize_kernal(coh);
        coh2d=coh;%
        
        if coh_dim == 3,coh=make_slice_3D(coh2d,zeros([kern_size,kern_size,kern_size]),orien0,(kern_size+1)/2);end
    
    end
    
    %if params.use_fftconv == 0,
    if ndims(coh) == 3,
        if kern_size >=0,coh=extract_max(coh,kern_size,kern_size,kern_size);else
            disp('Kern size = ROI size....'),end
    end
    %end
    
    if params.symmetrize_kernal == 1,
        disp('Symmetrizing coherence function....')
        coh=symmetrize_kernal(coh);
    end
    
    coh=abs(coh)/sum(sum(sum(abs(coh))));
    params.coh=coh;
    
end
%#########################################################################
%% Gauss kern minimizer.  works quite well    
if strcmp(params.pcdi_type,'gauss')
    
    disp('Updating coherence function using line search with gaussian model...')
    %if gauss_min == 1
    x0=params.kernalxy;
    
    %% when doing ,y,z do two runs of the minimisation
    %doing two iterations, this seems to be the most reliable.  one is fine
    %most of the time but run it twice should be good everytime
    if ndims(Mm) == 3
        sz=size(Mm);
        %do a 2d search first for the xy then use that as a gues for the z
        Mm_r0=zero_pad(squeeze(Mm(:,:,round(sz(3)/2))),params.roi(1),params.roi(2));
        Mk_r0=zero_pad(squeeze(Mk(:,:,round(sz(3)/2))),params.roi(1),params.roi(2));
        kernalxy = gauss_kern_minimizer_ver2(Mk_r0,Mm_r0,x0,params.angle,do(1),do(2),0,params.angle_not_do);
        %disp(kernalxy)
        kernalxy = gauss_kern_minimizer_ver2(Mk_r,Mm_r,kernalxy,params.angle,do(1),do(2),do(3),params.angle_not_do);
        %disp(kernalxy)
        kernalxy = gauss_kern_minimizer_ver2(Mk_r,Mm_r,kernalxy,params.angle,do(1),do(2),do(3),params.angle_not_do);
        coh=gauss_3D(gn,gn,gn,kernalxy(1),kernalxy(2),kernalxy(3));
        params.coh=coh;
    end
    %% do only once for 2D
    if ndims(Mm) ==2
        
        x0=determine_sigmas_in(orien0,x0);
        
        %disp(x0)
        [kernalxy angle]= gauss_kern_minimizer_ver2(Mk_r,Mm_r,x0,params.angle,do(1),do(2),do(3),params.angle_not_do);
        %[kernalxy angle]= gauss_kern_minimizer_ver2(Mk_r,Mm_r,kernalxy,angle,do(1),do(2),do(3),params.angle_not_do);
        
        disp(orien0)
        params.angle=angle;
        
        coh2d=gauss_2D(11,11,kernalxy(1),kernalxy(2),angle);
        
        %coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
        
        if coh_dim == 3,            %if data is 3d then check if 2d coh needs to be made 
            
            coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
            
        end
        
        
        if ndims(Mm) == 2, kernalxy=determine_sigmas_out(orien0,kernalxy);end
        
        params.coh=coh;
        
    end
   
    disp([kernalxy,params.angle])
    params.kernalxy=kernalxy;
        

end

%% simulated annealing for gaussian model
if strcmp(params.pcdi_type,'gauss_sa')
    disp('Updating coherence function using simulated annealing and a gaussian model...')
    x0=params.kernalxy;
    x0=determine_sigmas_in(orien0,x0);
    x0=[x0(1),x0(2)];
    
    if params.angle_not_do == 1,
        
        f=@(x)objfun_gauss2D(x,Mm_r,Mk_r,params.angle);
        x=simulannealbnd(f,x0,lb(1:2),ub(1:2));
    else
        x0=[x0(1),x0(2),params.angle];
        f=@(x)objfun_gauss2D_ang(x,Mm_r,Mk_r);
        x=simulannealbnd(f,x0,[lb(1:2),0],[ub(1:2),45]);
        params.angle=x(3);
    end
    
   
    if numel(x) == 2,coh2d=gauss_2D(11,11,x(1),x(2),parmas.angle);end
    if numel(x) == 3,coh2d=gauss_2D(11,11,x(1),x(2),x(3));end
    
    coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    
    params.coh=coh;
    
    kernalxy=determine_sigmas_out(orien0,[x(1),x(2)]);
    params.kernalxy=kernalxy;
    disp([kernalxy,params.angle])
    
    
    
end

%% pattern search for gaussian model.  Works quite well
if strcmp(params.pcdi_type,'gauss_ps')
    disp('Updating coherence function using pattern search and a gaussian model...')
    x0=params.kernalxy;
    options=psoptimset('Display','off','TolFun',1e-7,'TolMesh',1e-6,'MeshContraction',.5,'PollingOrder','Random','MeshExpansion',3);

    if ndims(Mm) == 2
    
        x0=determine_sigmas_in(orien0,x0);
        x0=[x0(1),x0(2)];
   
        if sum(params.do) == 2, %doing 1d
            if params.angle_not_do == 1,
                x0=[x0(1),x0(2)];

                if do(1) == 0,
                    f=@(x)objfun_gauss2D_x(x,Mm_r,Mk_r,params.angle);
                    x=patternsearch(f,x0(1),[],[],[],[],[lb(1)],[ub(1)],[],options);
                    x=[x(1),0];
                end
                if do(2) == 0,
                    f=@(x)objfun_gauss2D_y(x,Mm_r,Mk_r,params.angle);
                    x=patternsearch(f,x0(2),[],[],[],[],[lb(2)],[ub(2)],[],options);
                    x=[0,x(1)];
                end

            else
                x0=[x0(1),x0(2)];
                other=0.15;         %this is required to get the roated gaussian
                                    %if it =0, then doesn't work
                if do(1) == 0,
                    f=@(x)objfun_gauss2D_x_ang(x,Mm_r,Mk_r,other);
                    x=patternsearch(f,[x0(1),params.angle],[],[],[],[],[lb(1),-45],[ub(1),180],[],options);
                    params.angle=x(2);
                    x=[x(1),0];
                end
                if do(2) == 0,
                    f=@(x)objfun_gauss2D_y_ang(x,Mm_r,Mk_r,other);
                    x=patternsearch(f,[x0(2),params.angle],[],[],[],[],[lb(2),-45],[ub(2),180],[],options);
                    params.angle=x(2);
                    x=[0,x(1)];
                end

            end
        else   %assuming 2d and both paramters and maybe angle
            if params.angle_not_do == 1,

                f=@(x)objfun_gauss2D(x,Mm_r,Mk_r,params.angle);
                x=patternsearch(f,x0,[],[],[],[],[lb(1:2)],[ub(1:2)],[],options);
            else

                x0=[x0(1),x0(2),params.angle];
                f=@(x)objfun_gauss2D_ang(x,Mm_r,Mk_r);
                x=patternsearch(f,x0,[],[],[],[],[lb(1:2),0],[ub(1:2),180],[],options);
                params.angle=x(3);
            end
        end

        if numel(x) == 2,coh2d=gauss_2D(11,11,x(1),x(2),params.angle);end
        if numel(x) == 3,coh2d=gauss_2D(11,11,x(1),x(2),x(3));end

        %coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
        coh=coh2d;
        kernalxy=determine_sigmas_out(orien0,[x(1),x(2)]);
    end
    if ndims(Mm) == 3
        
        f=@(x)objfun_gauss3d(x,Mm_r,Mk_r); %was 3D
        x=patternsearch(f,x0,[],[],[],[],[lb],[ub],[],options);
        coh=gauss_3D(gn,gn,gn,x(1),x(2),x(3));
        kernalxy=x;
    end
    
    if coh_dim == 3,            %if data is 3d then check if 2d coh needs to be made 
        if sum(params.do(:)) ~= 0
            coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
        end
    end
    
    params.coh=coh;
    
    
    params.kernalxy=kernalxy;
    disp([kernalxy,params.angle])
    
    
    
end
%% pattern search for gaussian model.  Works quite well
if strcmp(params.pcdi_type,'kern_ps')
    disp('Updating coherence function using pattern search and a gaussian model with additional refinement....')
    x0=params.kernalxy;
    x0=determine_sigmas_in(orien0,x0);
    x0=[x0(1),x0(2)];
    
    options=psoptimset('Display','off','TolFun',1e-7,'TolMesh',1e-6,'MeshContraction',.5,'PollingOrder','Random','MeshExpansion',3);
   
    if sum(params.do) == 2, %doing 1d
        if params.angle_not_do == 1,
            x0=[x0(1),x0(2)];
            
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x(x,Mm_r,Mk_r,params.angle);
                x=patternsearch(f,x0(1),[],[],[],[],[lb(1)],[ub(1)],[],options);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y(x,Mm_r,Mk_r,params.angle);
                x=patternsearch(f,x0(2),[],[],[],[],[lb(2)],[ub(2)],[],options);
                x=[0,x(1)];
            end
           
        else
            x0=[x0(1),x0(2)];
            other=0.15;         %this is required to get the roated gaussian
                                %if it =0, then doesn't work
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x_ang(x,Mm_r,Mk_r,other);
                x=patternsearch(f,[x0(1),params.angle],[],[],[],[],[lb(1),-45],[ub(1),90],[],options);
                params.angle=x(2);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y_ang(x,Mm_r,Mk_r,other);
                x=patternsearch(f,[x0(2),params.angle],[],[],[],[],[lb(2),-45],[ub(2),90],[],options);
                params.angle=x(2);
                x=[0,x(1)];
            end
            
        end
    else   %assuming 2d and both paramters and maybe angle
        if params.angle_not_do == 1,

            f=@(x)objfun_gauss2D(x,Mm_r,Mk_r,params.angle);
            x=patternsearch(f,x0,[],[],[],[],[lb(1:2)],[ub(1:2)],[],options);
        else

            x0=[x0(1),x0(2),params.angle];
            f=@(x)objfun_gauss2D_ang(x,Mm_r,Mk_r);
            x=patternsearch(f,x0,[],[],[],[],[lb(1:2),0],[ub(1:2),90],[],options);
            params.angle=x(3);
        end
    end

    coh2d=gauss_2D(7,7,x(1),x(2),params.angle);
    %now update the values using a direct minimisation scheme
    options=psoptimset('Display','off','TolFun',1e-7,'TolMesh',1e-6,'MeshContraction',.5,'PollingOrder','Random','MeshExpansion',3);
    x0=coh2d;
    lb=zeros(size(x0));
    ub=ones(size(x0));
    f=@(x)objfun_kern2D(x,Mm_r,Mk_r);
    x=patternsearch(f,x0,[],[],[],[],lb,ub,[],options);
    coh2d=padarray(x,[2,2]);
    
    coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    
    params.coh=coh;
    
    kernalxy=determine_sigmas_out(orien0,[x(1),x(2)]);
    params.kernalxy=kernalxy;
    disp([kernalxy,params.angle])
    
    
    
end
%% unconstrained minimisation.  not as good as the others.
if strcmp(params.pcdi_type,'gauss_um')
    disp('Updating coherence function using unconstrained minimization and a gaussian model...')
    x0=params.kernalxy;
    x0=determine_sigmas_in(orien0,x0);

    options = optimset('Display','off');
    if sum(params.do) == 2, %doing 1d
        if params.angle_not_do == 1,
            x0=[x0(1),x0(2)];
            
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x(x,Mm_r,Mk_r,params.angle);
                x=fminunc(f,x0(1),options);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y(x,Mm_r,Mk_r,params.angle);
                x=fminunc(f,x0(2),options);
                x=[0,x(1)];
            end
           
        else
            x0=[x0(1),x0(2)];
            other=0.15;         %this is required to get the roated gaussian
                                %if it =0, then doesn't work
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x_ang(x,Mm_r,Mk_r,other);
                x=fminunc(f,[x0(1),params.angle],options);
                params.angle=x(2);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y_ang(x,Mm_r,Mk_r,other);
                x=fminunc(f,[x0(2),params.angle],options);
                params.angle=x(2);
                x=[0,x(1)];
            end
            
        end
    else   %assuming 2d and both paramters and maybe angle
        if params.angle_not_do == 1,

            f=@(x)objfun_gauss2D(x,Mm_r,Mk_r,params.angle);
            x=fminunc(f,x0,options);
        else

            x0=[x0(1),x0(2),params.angle];
            f=@(x)objfun_gauss2D_ang(x,Mm_r,Mk_r);
            x=fminunc(f,x0,options);
            params.angle=x(3);
        end
    end
   
    if numel(x) == 2,coh2d=gauss_2D(11,11,x(1),x(2),params.angle);end
    if numel(x) == 3,coh2d=gauss_2D(11,11,x(1),x(2),x(3));end
    
    coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    params.coh=coh;
    kernalxy=determine_sigmas_out(orien0,[x(1),x(2)]);
    params.kernalxy=kernalxy;
    disp([kernalxy,params.angle])
    
end
%% constrained minimisation
if strcmp(params.pcdi_type,'gauss_cm')
    disp('Updating coherence function using constrained minimization and a gaussian model...')
    x0=params.kernalxy;
    x0=determine_sigmas_in(orien0,x0);

    %3 algorithms to choose from
    %1. active-set -  not so good
    %2. interior-point - works quite well
    %3. sqp - not so good
    
    options = optimset('Display','off','Algorithm','interior-point');
    
     if sum(params.do) == 2, %doing 1d
        if params.angle_not_do == 1,
            x0=[x0(1),x0(2)];
            
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x(x,Mm_r,Mk_r,params.angle);
                x=fmincon(f,x0(1),[],[],[],[],[lb(1)],[ub(1)],[],options);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y(x,Mm_r,Mk_r,params.angle);
                x=fmincon(f,x0(2),[],[],[],[],[lb(2)],[ub(2)],[],options);
                x=[0,x(1)];
            end
           
        else
            x0=[x0(1),x0(2)];
            other=0.15;         %this is required to get the roated gaussian
                                %if it =0, then doesn't work
            if do(1) == 0,
                f=@(x)objfun_gauss2D_x_ang(x,Mm_r,Mk_r,other);
                x=fmincon(f,[x0(1),params.angle],[],[],[],[],[lb(1),-45],[ub(1),45],[],options);
                params.angle=x(2);
                x=[x(1),0];
            end
            if do(2) == 0,
                f=@(x)objfun_gauss2D_y_ang(x,Mm_r,Mk_r,other);
                x=fmincon(f,[x0(2),params.angle],[],[],[],[],[lb(2),-45],[ub(2),45],[],options);
                params.angle=x(2);
                x=[0,x(1)];
            end
            
        end
    else   %assuming 2d and both paramters and maybe angle
        if params.angle_not_do == 1,

            f=@(x)objfun_gauss2D(x,Mm_r,Mk_r,params.angle);
            x=fmincon(f,x0,[],[],[],[],[lb(1:2)],[ub(1:2)],[],options);
        else

            x0=[x0(1),x0(2),params.angle];
            f=@(x)objfun_gauss2D_ang(x,Mm_r,Mk_r);
            x=fmincon(f,x0,[],[],[],[],[lb(1:2),0],[ub(1:2),45],[],options);
            params.angle=x(3);
        end
     
    end
    
    if numel(x) == 2,coh2d=gauss_2D(11,11,x(1),x(2),params.angle);end
    if numel(x) == 3,coh2d=gauss_2D(11,11,x(1),x(2),x(3));end
    
    coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    params.coh=coh;
    kernalxy=determine_sigmas_out(orien0,[x(1),x(2)]);
    params.kernalxy=kernalxy;
    disp([kernalxy,params.angle])
   
end
%% least squares minimisation.  can do arbitary coh functions. not really sure how this performs
if strcmp(params.pcdi_type,'lsqs')
    disp('Updating coherence function using least-square minimization...')
    
    x0=extract_3D_slice(params.coh,orien0,6 );
    options = optimset('TolFun',1e-10,'Algorithm','trust-region-reflective','TolX',1e-8);
    lb=zeros(size(x0));
    ub=ones(size(x0));
    
    f=@(x)least_squares(x,Mm_r,Mk_r);
    [x resnorm residual]=lsqnonlin(f,x0,lb,ub,options);
   
    x=abs(x)/sum(sum(abs(x)));
    
    coh2d=x;
    coh2d=symmetrize_kernal(coh2d);
    coh=make_slice_3D(coh2d,zeros([11,11,11]),orien0,6);
    
    params.coh=coh;
    
end
%%
end

function [xyz]=determine_sigmas_out(orien0,kernalxy)
% determine what the order of the sigmas

if strcmp(orien0,'xy'), xyz=[kernalxy(1),kernalxy(2),0];end
if strcmp(orien0,'zx'), xyz=[kernalxy(2),0,kernalxy(1)];end
if strcmp(orien0,'xz'), xyz=[kernalxy(1),0,kernalxy(2)];end
if strcmp(orien0,'zy'), xyz=[0,kernalxy(2),kernalxy(1)];end


end

function [xyz]=determine_sigmas_in(orien0,kernalxy)
% determine what the order of the sigmas

if strcmp(orien0,'xy'), xyz=[kernalxy(1),kernalxy(2),0];end
if strcmp(orien0,'zx'), xyz=[kernalxy(3),kernalxy(1),0];end
if strcmp(orien0,'xz'), xyz=[kernalxy(1),kernalxy(3),0];end
if strcmp(orien0,'zy'), xyz=[kernalxy(3),kernalxy(2),0];end


end

function [coh3d] = make2Dcoh_3D(coh,orien0)

sz=size(coh);

newax=max(sz);

while mod(newax,2) == 0,newax=newax+1;end %make odd to embed in the center

cent=(newax+1)/2;%  get central pixel value of new axis

if strcmp(orien0,'xy'), 
    
    coh3d=zeros([sz,max(sz)]);
    coh3d(:,:,cent)=coh;
end    
    
if strcmp(orien0,'zx'), 
    %yxz
    %yx
    
    coh3d=zeros([max(sz),sz(1),sz(2)]);
    coh3d(cent,:,:)=coh;
end 
    
if strcmp(orien0,'xz'), 
    %yxz
    %yx
    
    coh3d=zeros([max(sz),sz(2),sz(1)]);
    coh3d(cent,:,:)=rot90(coh);
end 


if strcmp(orien0,'zy'), 
    %yxz
    %yx
    
    coh3d=zeros([sz(1),max(sz),sz(2)]);
   coh3d(:,cent,:)=coh;
end 



end

function params = set_determine_coh_defualts(params)

try
    params.min_sigmas;
catch
    params.min_sigmas=[0,0,0];
end
try
    params.max_sigmas;
catch
    params.max_sigmas=[2,2,2];
end

try
    params.lucy_its;
catch
    params.lucy_its=500;
end

try
    params.use_previous_coh;
catch
    params.use_previous_coh = 0;
end
try
    params.use_fftconv;
catch
    params.use_fftconv=0;
end
try
    params.pcdi_seperable;
catch
    params.pcdi_seperable=0;
end
try
    params.symmetrize_data;
catch
    params.symmetrize_data=0;
end
try
    params.symmetrize_kernal;
catch
    params.symmetrize_kernal=0;
end

end