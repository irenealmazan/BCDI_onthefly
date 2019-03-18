%% on the fly for loading tiff stack, remove aliens, put into data and run it quickly
clear all; close all; clc; %initialize_plot_PRL;

% load in tiff stack

cd ~/Documents/Hoydoo216bdsets/
curdir = pwd; %start on desktop

scannum2=[2623:2650];
viewdat =1;
cnt=1;
wantin = 1;
remalienfromprev = 1;
%
for   i=1:length(scannum2)
    close all;
    scannum = scannum2(i);
    
    %part a
    datdir = sprintf('~/Documents/Hoydoo216/Hoydoo216b_S%04d/',scannum);
    
    
    %part b
    %datdir = sprintf('~/Documents/Hoydoo216/Hoydoo216b_S%04d/',scannum);
    
    
		
try 
    cd(datdir);
    clear listOftifs;
		clear dir;
    listOftifs = dir('*.tif');
    nfilesmax = numel(listOftifs)
    ind = strcmp('ad_align_00000.tif',listOftifs(end).name);
    if ind ==1
        nfilesmax = nfilesmax-1;
    end
    
    ind = strcmp('ad_align_0000.tif',listOftifs(end).name);
    if ind ==1
        nfilesmax = nfilesmax-1;
    end
    
        ind = strcmp('ad_align_000.tif',listOftifs(end).name);
    if ind ==1
        nfilesmax = nfilesmax-1;
    end
    
if nfilesmax>25 && nfilesmax<150
    data = zeros(256,256,nfilesmax);
    cd(curdir);
    
    for nfiles = 0:nfilesmax-1
         I = imread([datdir listOftifs(nfiles+1).name]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir listOftifs(nfiles+1).name])
    end
    sname = sprintf('S%04d',scannum);
    
    if viewdat ==1
    play_data(data);
    end
    
    if wantin ==1
    resultal = input('require alien removal?');
    else
    resultal=0;
    end
   
    
if resultal == 1;
if remalienfromprev ==1
if cnt~=1
    for tcnt = 1:cnt-1
    data(abs(ycstore(tcnt)-10:ycstore(tcnt)+10)+1,abs(xcstore(tcnt)-10:xcstore(tcnt)+10)+1,abs(zcstore(tcnt)-2:zcstore(tcnt)+2))=0;
    end
end
end


		figure;

        
for gw=3:2:nfilesmax
			imagesc(log(data(:,:,gw)))
            title(['zsl ' num2str(gw) ' diff data for scannum' num2str(scannum)])
    %        imagesc(data(:,:,g))
			%pause(0.2)
            
            result = input('require alien removal?');
 
% 
                if result>=1
                    for cvg = 1:result
                disp('click on center of alien, surrounding 20 pixels will be zerod');
                [xc yc] = ginput(1);
                xc = round(xc);
                yc = round(yc);
                data(abs(yc-10:yc+10)+1,abs(xc-10:xc+10)+1,abs(gw-2:gw+2))=0;
                
                xcstore(cnt) = xc;
                ycstore(cnt) = yc;
                zcstore(cnt) = gw;
                
                cnt = cnt+1;
                
                clear xc yc;
                    end
                end

end

    clear gw;
    clf;
end
    
% if resultal == 1
%     play_data(data);
%    
%          close gcf;
% 		
% end 

 
 
   datasave = data(1:256,1:256,1:end);
    data = datasave;
    
    %for storing all the sets and then looking at their evolution
    %datast(:,:,:,cnt) = data;
    %cnt = cnt+1;
    

 %	 save(sname,'data');
 %disp('saved data')

%    pause(0.2);
%     if wantin ==1
%     result2 = input('appear phaseable');
%     else
%     result2 = 1;
%     end
   
%    if result2==1
 	 save(sname,'data');
     
%    else
        
%       scannum2(ind)=0; 
%    end


%     
%      result2 = input('appear phaseable');
%     if result2==1
%  	 save(sname,'data');
%     end
    
 
else ind= find(scannum==scannum2);
                        disp('nothing for')
                        scannum
			scannum2(ind)=0;

            
end
		cd(curdir);

catch me
    me
        disp('nothing for')
                        scannum
                        ind= find(scannum==scannum2);
			scannum2(ind)=0;
end
end
%clf;
scannum2(scannum2==0)=[];
%%
load(sname); %load it back up
%% pick some min data and zero it/pad it
data = sqrt(data);
mindata=3;
centdat=1;

ind=( data < mindata );
 data(ind)=0;
 
 data=fft_padjesse(data,[[1 1],1]);
 
 if centdat == 1
 data=center_array(data);
 disp('Centering data....')
 end


%% random guess
support_i=zero_pad_ver3(ones([round(.4*size(data,1)),round(size(data,2)*.4),round(size(data,3)*.4)]),size(data,1),size(data,2),size(data,3));
params.start_guess='random-data';        
sx=ceil(0.6*size(data,1));
sy=ceil(0.6*size(data,1));   
sz=ceil(0.6*size(data,3));   %x,y,z support size (pixels)

ss=round([sy,sx,sz]);            %its row column major so y is before x    
params.nn=size(data);


[rho support] = make_intial_support(data,ss,params);


%% shrinkwrap params

swsupport = 1;
thresh = 0.1;
sigthres = 1;

%% iteration time
Niter = 515;
algswitch = 40;
beta = .9; %HIO param
cnt = 1;
ERflag = 1;
%% iteration loop

for N=1:Niter %iteration loop
    
    if(~mod(N,algswitch)) ERflag = ERflag * -1;% display('switched alg'); 
    end
    
    
    %% modulus constraint
    Psi2 = fftshift(fftn((rho))); %guess of the obj at detector
    
    
    Psi_mod = data.* exp(1i*angle(Psi2)); %modulus projector
    

    Pmrho = (ifftn(fftshift(Psi_mod))); %back to real space
    
%    calculating the error

    %full calculation
    nume=sum(sum(sum(abs(abs(Psi2)-data).^2)));                 
    denom=sum(sum(sum(data.^2)));
    err(N)=nume/denom;   

     disp(['iteration number: (',num2str(N),'/',num2str(Niter),')','[',num2str(ERflag),']','  error =',num2str(err(N))])

    
    %try shrinkwrap here, see if it makes any difference
    if swsupport ==1
    if(~mod(N,10))
    %update support via shrinkwrap
    support=shrink_wrap(abs(Pmrho),thresh,sigthres,'gauss');
    
   [p1 p2 p3] = ind2sub(size(Pmrho(:,:,:,1)),find(abs(Pmrho(:,:,:,1))==max(max(max(abs(Pmrho(:,:,:,1)))))));
    p3=p3(end);
    
    %both
   subplot(3,1,1), semilogy(err);
    
    subplot(3,1,2), imagesc(squeeze(abs(Pmrho(:,p2(end),:,1)))); 
    subplot(3,1,3), imagesc(angle(Pmrho(:,:,p3,1)));
    drawnow;
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

[pnm yxz]=center_array(rho);
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
rho = rho/max(abs(rho(:)));

%%
figure; imagesc(abs(rho(:,:,end/2)));
figure; imagesc(angle(rho(:,:,end/2)));