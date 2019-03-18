%% do GA on the fly
% clear all; close all; clc;
% 
 addpath(genpath('/Users/ialmazn/Box Sync/forIrene/Matlab-Phasing-27May-2014'));
 addpath('/Users/ialmazn/Box Sync/forIrene/openspec-1.4/');
% addpath(genpath('/Users/ialmazn/Documents/MATLAB/ptycho/m_scripts/'))

% load('data_exp_struct.mat');
% 
% %%%% definition of detector experimental parameters
% pixsize = 55; %microns, for Merlin
% lam = etolambda(10400)*1e-4;
% Npix = 200;
% detdist = 0.529e6; %um - using scan#: 24915, 24918 w/ vert change of 0.5 deg. 
% d2_bragg = detdist * lam /(Npix*pixsize);
% depth = numel(data_exp);
% 
% %%%% geometrical parameters of the experiment
% th = 73.3;
% del = -32.6; %in plane
% gam = 0; %out of plane

%NW_diff_vectors_BCDI_v2;



%data = imgs;

% 
% %load('data_file2.mat');
% 
% %%{ This is specific from _v2
% 
% global  X Y Z d2_bragg ki_o kf_o
% 
% specnum = 92;
% 
% [specscan, errors] = openspec('Stephenson316a.spec',specnum );
% 
% % the geometry of the experiment:
% BCDI_diff_vectors;
% 
% 
% % load the tif files into a matrix:
% 
% filename =  ['Stephenson316a_S' num2str(specnum,'%04d') '_'];
% directname =  ['Stephenson316a_S' num2str(specnum,'%04d') '/'];
% 
% %[data_exp,rock,imgs] = BCDI_read_center_scans(filename,directname,specscan);
% [data_exp,rock,imgs] = BCDI_read_center_pad_scans(filename,directname,specscan);
% 
% %}

for ii = 1:numel(data_exp)    
    data(:,:,ii) = data_exp(ii).I;  
    data_exp(ii).dth_new = data_exp(ii).dth+data_exp(ii).dth_delta; 
    data_exp(ii).dqshift = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz]; 
 end

%data = imgs;

mindata=3;
centdat=0;
bindat = 0;

% ind=( data < mindata );
%  data(ind)=0;

 % padding the data matrix to 
% data=fft_padjesse(data,[[1 1],1]);
arrysize = 256;%200;
arrysize3=64;%51;

datap = zeros(arrysize,arrysize,arrysize3);

datap(arrysize/2-size(data,1)/2+1:arrysize/2+size(data,1)/2,arrysize/2-size(data,1)/2+1:arrysize/2+size(data,1)/2, ...
arrysize3/2-size(data,3)/2+1:arrysize3/2+size(data,3)/2) = data; %sample to the losless array size

data = datap;

 if centdat == 1
 data=center_array(data);
 disp('Centering data....')
 end
 
 if bindat ==1
     for qq=1:size(data,3)
        data_new(:,:,qq)=box_interp(data(:,:,qq),2,2,0);
     end
     
      data = data_new;
 end

 
 %data = sqrt(data);
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


[rho,support] = make_intial_support_Irene(data,ss,params);

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
    %rho_test = rho;
    
    %%%%% conventional fft
    %%{
    Psi2 = fftshift(fftn((rho))); %guess of the obj at detector
    
    
    Psi_mod = data.* exp(1i*angle(Psi2)); %modulus projector
    

    Pmrho = (ifftn(fftshift(Psi_mod))); %back to real space
    %}
    
    %%%%% slice ft
    
    %{
    [Psi2,Psi2_mod,rock_Psi2] = FT_Irene(data,data_exp,rho);
    
    Psi_mod = data.* exp(1i*angle(Psi2)); %modulus projector
    
    [Pmrho] = InvFT_Irene(data,data_exp,Psi_mod);
    
    %}
%    calculating the error
    if(~mod(N,20))
    %full calculation
    nume=sum(sum(sum(abs(abs(Psi2)-data).^2)));                 
    denom=sum(sum(sum(data.^2)));
    err(N)=nume/denom;   

     disp(['iteration number: (',num2str(N),'/',num2str(Niter),')','[',num2str(ERflag),']','  error =',num2str(err(N))])
    end
    %remove the twin
%     if N == twinit
%         support(:,:,:,qq) = support0.*support(:,:,:,qq);
%         disp('twin it')
%     end
    
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
figure; 
subplot(121);
imagesc(abs(rhobest(:,:,round(numel(data_exp)/2))));
%imagesc(abs(rho(:,:,round(numel(data_exp)/2))));
axis image;

subplot(122);
imagesc(angle(rhobest(:,:,round(numel(data_exp)/2))));
%imagesc(angle(rho(:,:,round(numel(data_exp)/2))));
axis image;

figure; 
subplot(121);
imagesc(abs(squeeze(rhobest(round(Npix/2),:,:))));
%imagesc(abs(squeeze(rho(round(Npix/2),:,:))));
axis image;

subplot(122);
imagesc(angle(squeeze(rhobest(round(Npix/2),:,:))));
%imagesc(angle(squeeze(rho(round(Npix/2),:,:))));
axis image;

if exist('NW','var')
    figure;
    subplot(121);
    imagesc(abs(squeeze(NW(round(Npix/2),:,:))));
    axis image;
    
    subplot(122);
    imagesc(angle(squeeze(NW(round(Npix/2),:,:))));
    axis image;
    
    figure; 
    subplot(121);
    imagesc(abs(NW(:,:,round(numel(data_exp)/2))));
    axis image;
    
    subplot(122);
    imagesc(angle(NW(:,:,round(numel(data_exp)/2))));
    axis image;
end
