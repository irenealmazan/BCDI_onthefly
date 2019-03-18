function [ prtf ] = prtf_MatlabPhasing(ddir)
%jclark
%ddir is for the reconstruction

nnames=rdirnames([ddir,'/*PARAMS.mat']); %load the params

load(nnames{1});

[ data ] =load_rec_data(params);


%pn=load_rec_from_dir(ddir,1);
nnames=rdirnames([ddir,'/*pn_g.mat']); %load the params
load(nnames{1});
[pn_g] = align_iterates(pn_g);
pn=sum(pn_g,4);

Ipnm=convn(abs(fftshift(fftn(fftshift(pn)))).^2,params.coh,'same');

Ipnm(Ipnm <0)=0;

%pp=0;
%for qq=1:2:size(Ipnm,3)
%    pp=pp+1;
%    Ipnmn(:,:,pp)=sum(Ipnm(:,:,qq:qq+1),3);
%    datan(:,:,pp)=sum(data(:,:,qq:qq+1),3);
%end
Ipnmn=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(Ipnm))),size(Ipnm,2),size(Ipnm,1),62)))));
datan=abs(fftshift(fftn(fftshift(zero_pad_ver3(fftshift(ifftn(fftshift(data))),size(Ipnm,2),size(Ipnm,1),62)))));

Ipnm=Ipnmn;
data=datan;

%get the 'shells'
nd=max(size(data));
Ipnm=zero_pad_ver3(Ipnm,nd,nd,nd);
data=zero_pad_ver3(data,nd,nd,nd);

[h k l]=register_3d_reconstruction(abs(Ipnm),abs(data));

Ipnm=circshift(Ipnm,-round([h,k,l]));

xyz=center_of_mass(data);
data=circshift(data,-round([xyz(2),xyz(1),xyz(3)]));
Ipnm=circshift(Ipnm,-round([xyz(2),xyz(1),xyz(3)]));



Ipnm=Ipnm./sum(Ipnm(:))*sum(data(:));

rad=1:(round(nd(1)));

Ipnm=sqrt(Ipnm);
data=sqrt(data);

for qq=2:numel(rad)
    
    [ circ ] = generate_sphere(nd(1),rad(qq))-generate_sphere(nd(1),rad(qq-1));

    circ=circ.*(data > 0);
    ind=find(circ == 1);
    
    %Intms(qq)=sum(data(ind));
    
    %Intks(qq)=sum(Ipnm(ind));
    Intms(qq)=mean(data(ind));
    
    Intks(qq)=mean(Ipnm(ind));
    
    
end

prtf=(Intks./Intms);

prtf(isnan(prtf))=0;
prtf(prtf == inf)=0;

%prtf=prtf./max(prtf(:));

prtf(prtf>1)=0;

end

