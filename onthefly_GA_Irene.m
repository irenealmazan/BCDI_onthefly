%% do GA on the fly
%clear all; close all; clc;

addpath(genpath('/Users/ialmazn/Box Sync/forIrene/Matlab-Phasing-27May-2014'));
addpath('/Users/ialmazn/Box Sync/forIrene/openspec-1.4/');
addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))

load('data_exp_struct.mat');

%data = imgs;


% update data structure & prepare initial guess:
 %%{
 for ii = 1:numel(data_exp)
     
    data(:,:,ii) = data_exp(ii).I;

    %inibeam(ii).A = probe;
    %inibeam(ii).A = circshift(probe, round([shift_guess(ii,2)  shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

    %%% THETA POSITIONS
    
    %uncorrect theta positions
    %%{
    dth_initials(ii) = data_exp(ii).dth;
    data_exp(ii).dth_new = data_exp(ii).dth; 
    data_exp(ii).dqshift = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz]; 
    %}
    
    %correct theta positions
    %{
     data_exp(ii).dqshift =  data_exp(ii).dqshift_delta; 
    %}
    
    
 end

 % Initialize some parameters:
 
 err = [];
 
mindata=3;
centdat=0;
bindat = 0;

ind=( data < mindata );
 data(ind)=0;



% set up experimental details
%thBragg2 = 34; % in deg, taken from labbook page 19
[Npix,Npiy,depth] = size(data);
%d2_bragg = camdist * params.lam/(Npix*params.det_px);

arrysize = Npix;
arrysize3 = depth;

[X,Y,Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
                    [-Npiy/2:Npiy/2-1]*d2_bragg,...
                    [-depth/2:depth/2-1]*d2_bragg);





 if centdat == 1 % done already in BCDI_read_center_pad_scans
 data=center_array(data);
 disp('Centering data....')
 end
 
 if bindat ==1
     for qq=1:size(data,3)
    data_new(:,:,qq)=box_interp(data(:,:,qq),2,2,0);
     end
     
      data = data_new;
 end

 
 %data = sqrt(data); % done already in BCDI_read_center_pad_scans
 
 figure(5);
 clf;
%% start the generation loop
for ng = 1:1

%% start the population loop
    for np=1:1


    %% random guess
     support1 = ones([round(.4*size(data,1)),round(size(data,2)*.4),round(size(data,3)*.4)]);
    [support_i] = zero_pad_ver_Irene(arrysize,arrysize,arrysize3, support1);

    %support_i=zero_pad_ver3(ones([round(.4*size(data,1)),round(size(data,2)*.4),round(size(data,3)*.4)]),size(data,1),size(data,2),size(data,3));
    params.start_guess='random-data';        
    sx=ceil(0.6*size(data,1));
    sy=ceil(0.6*size(data,1));   
    sz=ceil(0.6*size(data,3));   %x,y,z support size (pixels)

    ss=round([sy,sx,sz]);            %its row column major so y is before x    
    params.nn=size(data);


    [rho_ini,support] = make_intial_support_Irene(data,ss,params);

    rho = 1e-3.*rho_ini;
    %% breed in the "best one"

    if ng>1
       rho = sqrt(rho.*rhobest);
       disp('breeding in best sharpness')
    end


    %% shrinkwrap params

    swsupport = 1;
    thresh = 0.1;

    %% iteration time
    Niter = 515;
    algswitch = 40;


    beta = .9; %HIO param
    twinit = 45; %cut off half the array at this iteration
    checkforcnj = 90; %should be 90 or so
    checkagain = 290;

    % support0 = zeros(size(support_i));
    % for i=1:size(data,3)
    % support0(:,:,i)=tril(support_i(:,:,i));            %for twin removal
    % end

    cnt = 1;
    ERflag = 1;
    %% iteration loop

        for N=1:Niter %iteration loop

            if(~mod(N,algswitch)) ERflag = ERflag * -1;% display('switched alg'); 
            end


            %% modulus constraint
            Psi2 = fftshift(fftn((rho))); %guess of the obj at detector

            %[Psi2,Psi2_mod,rock_Psi2] = FT_Irene(data,data_exp,rho);
            
            Psi_mod = data.* exp(1i*angle(Psi2)); %modulus projector

            rock_Psi_mod = sum(sum(sqrt(Psi_mod.*conj(Psi_mod)),1),2);

            Pmrho = (ifftn(fftshift(Psi_mod))); %back to real space

           %[Pmrho] = InvFT_Irene(data,data_exp,Psi_mod);

            for jj=1:numel(data_exp)               
                   subplot(1,5,3); imagesc(sqrt(data_exp(jj).I)); axis image;title('Experimental dp');colorbar;
                   %subplot(1,5,4); imagesc(sqrt(Psi_mod(:,:,jj).*conj(Psi_mod(:,:,jj))));axis image;title(['Calculated dp' num2str(jj)]);
                   subplot(1,5,4); imagesc(sqrt(Psi2(:,:,jj).*conj(Psi2(:,:,jj))));axis image;title(['Calculated dp' num2str(jj)]);colorbar;
                   thscantoplot(jj) = data_exp(jj).dth;
                   %subplot(1,5,5);plot(thscantoplot(1:jj),rock(1:jj),'-ob');hold on; plot(thscantoplot(1:jj),rock_Psi2(1:jj),'xr');
                   subplot(1,5,5); plot(thscantoplot(1:jj),rock(1:jj),'-ob');hold on; plot(thscantoplot(1:jj),squeeze(rock_Psi2(1:jj)),'xr')
                   
                   drawnow;
            end
            
             
           %try shrinkwrap here, see if it makes any difference
           if swsupport ==1
               if(~mod(N,5))
                   %update support via shrinkwrap
                   support=shrink_wrap(abs(Pmrho),thresh,1,'gauss');

               end
           end

           %% apply constraint in real space
           if ERflag ==1

               rho = support.*Pmrho;
           else
               %HIO:
               rho = support.*Pmrho + (1-support).*(rho-beta*Pmrho);

           end

            %calculating the error
            if(~mod(N,1))
                %full calculation
                nume=sum(sum(sum(abs(abs(Psi2)-data).^2)));                 
                denom=sum(sum(sum(data.^2)));
                err(N)=nume/denom;   

                 disp(['iteration number: (',num2str(N),'/',num2str(Niter),')','[',num2str(ERflag),']','  error =',num2str(err(N))])
            end

            size_rho_half = round(size(rho,3)/2);
            subplot(1,5,1); imagecomp(rho(:,:,size_rho_half)); colorbar; axis image; %zoom(1.5);
            subplot(1,5,2); plot(log10(err));
           
            drawnow;
        %   
            %remove the twin
        %     if N == twinit
        %         support(:,:,:,qq) = support0.*support(:,:,:,qq);
        %         disp('twin it')
        %     end

           


        end %end of iteration loop

    [pnm yxz]=center_array((rho));
    support=circshift(support,yxz);

    %do COM centering, will already be close
    xyz=center_of_mass(abs(pnm).*support);
    xyz=-1*round([xyz(2),xyz(1),xyz(3)]);
    pn=circshift(real(pnm),xyz)+1i*circshift(imag(pnm),xyz);
    support=circshift(support,xyz);

    disp('removing phase offset....')
    disp('')
    sz=size(pn);
    i=round(sz(1)/2);
    j=round(sz(2)/2);
    k=round(sz(3)/2);
    phi0=atan2(imag(pn(i,j,k)),real(pn(i,j,k)));
    pn=pn*exp(-1i*phi0);

    rho=remove_ramp_pn_ups(pn,3);

    rhostore(:,:,:,np) = rho;
    errstor(np) = err(end);
    sharpstore(np) = sum(abs(rho(:)).^4); %want min of this

    clear rho support;
    end

ind = find(min(sharpstore)==sharpstore);

rhobest = rhostore(:,:,:,ind);

disp('find new rhobest')

%% then want to choose minimum, use as seed for next reconstruction
%rho_n' = sqrt(rho_n * rho_alpha) rho_n = random again, rho_alpha = best
%sharpness

end

disp('final chi metric')
errstor(ind)
%%
rhobest = rhobest/max(abs(rhobest(:)));
figure; imagesc(abs(rhobest(:,:,end/2)));
figure; imagesc(angle(rhobest(:,:,end/2)));