%generating comparison scans
%load in the lab.mat file

%%
clc; clear all;

sname = 40;
com2 = 'A';

numscans=5;
iterates = zeros(256,256,48,numscans);
cohers = zeros(size(iterates));
linup = 1; %try to shift to match


figure
hold on
%do the loading into a 4D matrix
for i=1:numscans
scannum=i;

datdir = ['/Users/a4894z/Dropbox/Matlab_phasing_AUG13/GA/' int2str(sname) com2 sprintf('r%01d',scannum) '/'];

%fname = [sprintf('%02dAr%01d',sname,scannum) '-LAB.mat']; 

fname = [sprintf('%02d',sname) com2 sprintf('r%01d',scannum) '-LAB.mat']; 

fname2 = [sprintf('%02d',sname) com2 sprintf('r%01d',scannum) '-LAB-C.mat'];


anew = load([datdir fname]); %load amp/phase
cnew= load([datdir fname2]); %load coh

iterates(:,:,:,i) = anew.array;

cohers(:,:,:,i) = cnew.array;

subplot(5,2,i); imagesc(abs(anew.array(:,:,48/2))); zoom(3)
subplot(5,2,i+numscans); imagesc(abs(cnew.array(:,:,48/2))); zoom(3)


end

[ylim xlim zlim] = size(anew.array);


%now do the averaging
ref = iterates(:,:,:,1);

for qq=2:numscans      
        
   
    cnj_rf=is_conj_ref(ref,iterates(:,:,:,qq))
            
    if cnj_rf ==1
        iterates(:,:,:,qq)=conj_reflect(iterates(:,:,:,qq));
    end 
    
    %try to line them up using essentially cross correlation
    if linup == 1
    [h k l] = register_3d_reconstructions(abs(ref),abs(iterates(:,:,:,qq)));
    next_c = sub_pixel_shift(iterates(:,:,:,qq),h,k,l);
    iterates(:,:,:,qq) = next_c;
    end
    
end

asum = mean(iterates,4); %average amp and phase at same time

%save the average in a mat file
%v43 = genvarname([int2str(sname) com2]);
%eval([v43 '= asum;']);
%save([int2str(sname) com2],v43); %saves it as 40A.mat, x40A is asum


%%
%assume you loaded all the arrays, saved like x40A

subplot(3,2,1); imagesc(abs(x40A(:,:,48/2)).*(abs(x40A(:,:,48/2))>=0.2)); zoom(2)
subplot(3,2,2); imagesc(angle(x40A(:,:,48/2))); zoom(2)

subplot(3,2,3); imagesc(abs(x40B(:,:,48/2)).*(abs(x40B(:,:,48/2))>=0.2)); zoom(2)
subplot(3,2,4); imagesc(angle(x40B(:,:,48/2))); zoom(2)

subplot(3,2,5); imagesc(abs(x40C(:,:,48/2)).*(abs(x40C(:,:,48/2))>=0.2)); zoom(2)
subplot(3,2,6); imagesc(angle(x40C(:,:,48/2))); zoom(2)

subplot(9,2,7); imagesc(abs(x40D(:,:,48/2))); zoom(2)
subplot(9,2,8); imagesc(angle(x40D(:,:,48/2))); zoom(2)

subplot(9,2,9); imagesc(abs(x40E(:,:,48/2))); zoom(2)
subplot(9,2,10); imagesc(angle(x40E(:,:,48/2))); zoom(2)

subplot(9,2,11); imagesc(abs(x40F(:,:,48/2))); zoom(2)
subplot(9,2,12); imagesc(angle(x40F(:,:,48/2))); zoom(2)

subplot(9,2,13); imagesc(abs(x40G(:,:,48/2))); zoom(2)
subplot(9,2,14); imagesc(angle(x40G(:,:,48/2))); zoom(2)

subplot(9,2,15); imagesc(abs(x40H(:,:,48/2))); zoom(2)
subplot(9,2,16); imagesc(angle(x40H(:,:,48/2))); zoom(2)

subplot(9,2,17); imagesc(abs(x40I(:,:,48/2))); zoom(2)
subplot(9,2,18); imagesc(angle(x40I(:,:,48/2))); zoom(2)









%%
%looks like good cutoff might be 0.25
%check and see what the cutoff does to the amplitude
pval=0.2;
figure;
for zs=1:zlim
imagesc(abs(asum(:,:,zs).*(abs(asum(:,:,zs))>=pval)))
pause
end

pval=0.2;
figure;
for zs=1:zlim
imagesc(angle(asum(:,:,zs).*(abs(asum(:,:,zs))>=pval)))
pause
end

xmax = xlim;
ymax = ylim;
zmax = zlim;

pval = 0.2; %stick with this for now
phsf = 2*pi*sqrt(3)/8.17; %in angstrom

phasesp = (angle(asum).*(abs(asum)>=pval))/phsf; 
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,180,180,350); %use spacing of 17nm
clear uy uz;
uxnans = ux;
cdata = smooth3(ux,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
cdata(isnan(cdata)==1) = 0;
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);



%visualize the strain
figure;
for zsl = zlim/2
imagesc(cdata(:,:,zsl)*10^3); %multiply to make the scalebar easy to read
%zoom(4)
%caxis([-.25 .25])
%pause
end

figure;
%visualize the phase
for zsl = 1:zlim
imagesc(phasesp(:,:,zsl));
%zoom(4)
%caxis([-1 1])
pause
end

%do the isosurface projection, do before you run initialize plot

pval = 0.03; %stick with this for now
phsf = 2*pi*sqrt(3)/8.17; %in angstrom

phasesp = (angle(asum).*(abs(asum)>=pval))/phsf; 
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,180,180,350); %use spacing of 17nm
clear uy uz;
uxnans = ux;
cdata = smooth3(ux,'box',[3 3 1]);
%cdata = smooth3(phasesp,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
cdata(isnan(cdata)==1) = 0;
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);


figure;
p = patch(isosurface(x,y,z,abs(asum),0.2));
isonormals(x,y,z,abs(asum),p);
isocolors(x,y,z,cdata*10^3,p);
set(p,'FaceColor','interp','EdgeColor','none')
%daspect([1 1 1]);
axis off
camlight(230,250); 
lighting GOURAUD; 
%caxis([-1 1])
title('402 strain proj')
view(3); 



%pictures are good but numbers are better. Break up like before
%and calc total/surf/core rms strain
pval = 0.2;
phsf = 2*pi*sqrt(3)/8.17; %in angstrom
[ylim xlim zlim] = size(asum);
xmax = xlim;
ymax = ylim;
zmax = zlim;

phasesp = (angle(asum).*(abs(asum)>=pval))/phsf; 
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,170,170,350); %use spacing of 17nm
clear uy uz;
uxnans = ux;
cdata = smooth3(ux,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
cdata(isnan(cdata)==1) = 0;
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);

%center
uxcent = (ux.*(abs(asum)>=0.6));
uxcent(uxcent==0)=nan;
%imagesc(uxcent(:,:,zlim/2));


%edge is between .2 and 0.25
uxedge = (ux.*(abs(asum)<=0.45));
uxedge(uxedge==0)=nan;
%imagesc(uxedge(:,:,zlim/2));


%calculate the quant metrics
uxqtot = reshape(ux,xmax*ymax*zmax,1);
uxqcent = reshape(uxcent,xmax*ymax*zmax,1);
uxqedge = reshape(uxedge,xmax*ymax*zmax,1);


uxnmean = nanmean(uxqtot);
uxnmeancent = nanmean(uxqcent);
uxnmeanedge = nanmean(uxqedge);

%average of abs
uxavgtot = nanmean(abs(uxqtot));
uxavgcent = nanmean(abs(uxqcent));
uxavgedge = nanmean(abs(uxqedge));
%standard deviation
uxstdtot = nanstd(uxqtot);
uxstdcent = nanstd(uxqcent);
uxstdedge = nanstd(uxqedge);
%maximums
uxmaxtot = max(abs(uxqtot)); %this also might be useful
uxmaxcent = max(abs(uxqcent));
uxmaxedge = max(abs(uxqedge));
%rms
uxrmstot = sqrt(nanmean(uxqtot.^2));
uxrmscent = sqrt(nanmean(uxqcent.^2));
uxrmsedge = sqrt(nanmean(uxqedge.^2));

save('40Bavgd','asum','ssum','csum','uxrmstot','uxrmscent','uxrmsedge')


[xmax ymax zmax] = size(asum);
figure;
subplot(2,2,1); imagesc(abs(asum(:,:,zmax/2))); title('Amplitude center slice');
zoom(3);

subplot(2,2,2); imagesc(angle(asum(:,:,zmax/2))); title('Phase center slice');
zoom(3);

subplot(2,2,3); imagesc(cdata(:,:,zmax/2)); title('Strain center slice');
zoom(3); colorbar;
pval = 0.03; %stick with this for now
phsf = 2*pi*sqrt(3)/8.17; %in angstrom

phasesp = (angle(asum).*(abs(asum)>=pval))/phsf; 
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,180,180,350); %use spacing of 17nm
clear uy uz;
uxnans = ux;
cdata = smooth3(ux,'box',[3 3 1]);
%cdata = smooth3(phasesp,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
cdata(isnan(cdata)==1) = 0;
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
subplot(2,2,4); p = patch(isosurface(x,y,z,abs(asum),0.25));
isonormals(x,y,z,abs(asum),p);
isocolors(x,y,z,cdata*10^3,p);
set(p,'FaceColor','interp','EdgeColor','none')
%daspect([1 1 1]);
axis off
camlight(230,250); 
lighting GOURAUD; 
%caxis([-1 1])
title('Strain proj')

