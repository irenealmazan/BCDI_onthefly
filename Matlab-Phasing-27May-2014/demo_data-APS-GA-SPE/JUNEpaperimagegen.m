%June paper images. Keep up to date with which parameters are used for what
%use other code to test ideas, not this one until you're sure it's better
exp4ydat = [4.7,4.6, 4.52, 4.25, 4.06, 3.57, 4.06, 4.47, 4.67,4.74, 4.77, 4.85...%first curve
    3.95,4.735,3.46,4.72,3.98...%cycle a bit, then run cycles 6-10 w no scans, then
   3.57, 4.08, 4.42, 4.68, 4.74, 4.68, 4.6, 4.42, 3.95,3.66...%final curve
   ];

intpts = [4.7 3.5 4.7 3.5 4.7 3.5 4.7 3.5 4.7];

zpts = zeros(size(intpts)); %add this to the things I calculate

capdat = [1 6 12 13 14 15 16 17 18 22 27]; %list of full charge/full discharge states

exp4ddat=[
0.53844
0.21515
0.1786
0.14955
0.11008
1.39e-04
0.05268
0.11342
0.37208
0.72746
0.86812
0.89536
0.06519
0.55164
0.03544
0.50202
0.07728
1.39e-04
0.05893
0.10702
0.37208
0.72746
0.49925
0.21515
0.16484
0.06519
0.03878];

%use this for the lattice const graphs, stick with aeff for actual strain
%calcs
aeff4=[
    8.098841292
8.082505276
8.091335852
8.144877527
8.137046808
8.183136774
8.190092491
8.177716293
8.177138676
8.142470778
8.155375217
8.142351682
8.140935712
8.141033592
8.187651256
8.172903312
8.132190681
8.173910569
8.180330397
8.175844893
8.169025791
8.114736231
8.098563875
8.143313028
8.174397507
8.193170026
8.200700122];

%say you want to add in strain from average lattice constant shift
%relative to the previous one
avgstrainrel=diff(aeff4)./aeff4(1:end-1);

avgstrainabs=(aeff4-aeff4(1))./aeff4(1); %relative to first point

%josh PXRD data, compare with pts 6:12 in my data
josh_delta=[
0
0.21307
0.41
0.49
0.64254
0.7384
0.841
0.888];

josh_lat1=[
	8.1782
8.177263
8.16835
8.1405
8.139617
8.146
8.132
8.138];

josh_lat2=[
	8.0904
8.089
8.0747
8.0783];

josh_lat3=[
	8.00195
8.000195
8.00538];

%josh PXRD data for discharge
josh_delta_d=[0.809, 0.548, 0.486, 0.361, 0.295, 0.163, .052];

josh_lat1_d=[8.120433, 8.1286, 8.1283, 8.1344, 8.13,8.156, 8.166];
    
josh_lat2_d=[8.0904,8.0939]; 

josh_lat3_d=[8.0127];
	
	

aeff=[
8.1275
    8.1272
    8.1454
    8.1885
    8.1884
    8.2126
    8.2127
    8.2036
    8.1945
    8.1736
    8.1738
    8.1644
    8.1901
    8.1738
    8.2277
    8.2076
    8.2017
    8.2372
    8.2374
    8.2374
    8.2284
    8.1969
    8.1883
    8.2320
    8.2650
    8.2754
    8.2879];

%for making surface strain propagation
%you will need exp4alldat and aeffexp4
avgamp = mean(ampexp4(:,:,:,1:14),4); %use first 12 bc all have 64 z dim
avgamp2 = mean(ampexp4(:,:,:,[15:17 19:20 22:27]),4);
avgamp3 = mean(ampexp4(:,:,:,[18 21]),4);
%1,6,12-14 use avgamp, 15 16 17 22 27 use avgamp2, 18 uses its own
szvec = [64 64 64 64 64 64 64 64 64 64 64 64 64 64 96 96 96 80 96 96 80 96 96 96 96 96 96]; 

%SURFACE STRAIN PROPAGATION
for q=[1:6]
pval2 = 0.03; %doesn't matter so much, use .03 or 0.1
%aeff2 = aeff(q); %this changes for every scan, put in Angstrom
aeff2 = aeff4(q);
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

if szvec(q) == 64
    ampp = avgamp;
elseif szvec(q) == 96
    ampp = avgamp2;
elseif szvec(q) == 80
    ampp = avgamp3;
end

phasesp = (phexp4(:,:,:,q).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
ux(isnan(ux)==1) = 0;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);

%ux = (ux + avgstrainabs(q))*100; %put the strain into percentage

figure;
cdata = smooth3(ux,'box',[3 3 1]);
p = patch(isosurface(x,y,z,ampp,0.2));
isonormals(x,y,z,ampp,p);
isocolors(x,y,z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none')
%daspect([1 1 1]);
axis tight
camlight(230,250); 
lighting GOURAUD; 
%title('4.7 1st Charge ','FontSize',16)
colorbar;
%set(gca, 'ydir','reverse')
view(3); 
axis off;
caxis([-5e-4 5e-4]); 
%caxis([-1 1]);
set(gca,'fontsize',16)
end


%CORE STRAIN PROPAGATION use 0.8 or 0.5 gotta decide
for q=[1 3]
pval2 = 0.1; %doesn't matter so much, use .03 or 0.1
aeff2 = aeff(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

%ampp = ampexp4(:,:,:,q);
ampp = avgamp;
phasesp = (phexp4(:,:,:,q).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,200,200,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
ux(isnan(ux)==1) = 0;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);

figure;
cdata = smooth3(ux,'box',[3 3 1]);
p = patch(isosurface(x,y,z,ampp,0.53));
%p = patch(isosurface(x,y,z,ampp,0.6));
isonormals(x,y,z,ampp,p);
isocolors(x,y,z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none')
%daspect([1 1 1]);axis tight
camlight(230,250); 
lighting GOURAUD; 
%title('27 Discharge ','FontSize',16)
%colorbar;
%set(gca, 'ydir','reverse')
view(3); 
axis off;
caxis([-5e-4 5e-4]); 
set(gca,'fontsize',16)
end

%%
%SLICES OF STRAIN
for q=[1:6]
pval2 = 0.1; %doesn't matter so much, use .03 or 0.1
aeff2 = aeff4(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

if szvec(q) == 64
    ampp = avgamp;
elseif szvec(q) == 96
    ampp = avgamp2;
elseif szvec(q) == 80
    ampp = avgamp3;
end

phasesp = (phexp4(:,:,:,q).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
uxnans = ux;

%ux(isnan(ux)==0) = (ux(isnan(ux)==0) + avgstrainrel(q))*100; %put in the average strain

cdata = smooth3(ux,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
cdata(isnan(cdata)==1) = 0;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);


% figure;
% imagesc(uxx(:,:,96/2));
% colorbar;
% %caxis([-5e-4 5e-4]) %used to use -5e-5 5e-4
% zoom(5)
% axis off;


figure;
isosurface(cdata,1e-4)
%imagesc(cdata(:,:,96/2));
%colorbar;
%caxis([-5e-4 5e-4]) %used to use -5e-5 5e-4
%zoom(5)
%axis off;


%bin data
uxnans(uxnans > 1e-4) = 1e-4; 
uxnans(uxnans < -1e-4) = 1e-4;
ux(ux > 1e-4) = 1e-4; 
ux(ux < -1e-4) = 1e-4;

% for j=48-5:1:48+5
% figure;
% imagesc(ux(:,:,48));
% caxis([-1e-4 1e-4]) %used to use -5e-5 5e-4
% zoom(4)
% axis off;
% pause
% end

uxnansedit = reshape(uxnans,xlim*ylim*zlim,1);

%it looks like this way will work
varplot(q) = nanstd(uxnansedit);

end

%%

%QUANT METRICS AS VOLTAGE
for q=1:27

clear phasesp ampp ux;
%strain map, was using 0.1 before
pval=0.2;
aeff2 = aeff(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

ampp = ampexp4(:,:,:,q);
[xmax ymax zmax] = size(ampp);
phasesp = (phexp4(:,:,:,q).*(ampp>=pval))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan; %only want gradients where you have data
[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom
clear uy uz;


%way that gives interpretable quant metrics
%center is between 0.53 and 1
uxcent = (ux.*(ampp>=0.53));
uxcent(uxcent==0)=nan;

%check number of elements in your calculation
uxcent2 = uxcent;
uxcent2(isnan(uxcent2)==1) = 0; nnz(uxcent2);
clear uxcent2;


%edge is between .2 and 0.25
uxedge = (ux.*(ampp<=0.25));
uxedge(uxedge==0)=nan;

uxedge2 = uxedge;
uxedge2(isnan(uxedge2)==1) = 0; nnz(uxedge2);
clear uxedge2;

%calculate the quant metrics
uxqtot(:,q) = reshape(ux,xmax*ymax*zmax,1);
uxqcent(:,q) = reshape(uxcent,xmax*ymax*zmax,1);
uxqedge(:,q) = reshape(uxedge,xmax*ymax*zmax,1);


uxnmean(q) = nanmean(uxqtot(:,q));
uxnmeancent(q) = nanmean(uxqcent(:,q));
uxnmeanedge(q) = nanmean(uxqedge(:,q));

%average of abs
uxavgtot(q) = nanmean(abs(uxqtot(:,q)));
uxavgcent(q) = nanmean(abs(uxqcent(:,q)));
uxavgedge(q) = nanmean(abs(uxqedge(:,q)));
%standard deviation
uxstdtot(q) = nanstd(uxqtot(:,q));
uxstdcent(q) = nanstd(uxqcent(:,q));
uxstdedge(q) = nanstd(uxqedge(:,q));
%maximums
uxmaxtot(q) = max(abs(uxqtot(:,q))); %this also might be useful
uxmaxcent(q) = max(abs(uxqcent(:,q)));
uxmaxedge(q) = max(abs(uxqedge(:,q)));
%rms
uxrmstot(q) = sqrt(nanmean(uxqtot(:,q).^2));
uxrmscent(q) = sqrt(nanmean(uxqcent(:,q).^2));
uxrmsedge(q) = sqrt(nanmean(uxqedge(:,q).^2));

q

end

uxrmstot(2) = 207/1000000;

figure;
hold on
plot(exp4ddat(1:12),uxrmsedge(1:12),'b.')
plot(exp4ddat(18:27),uxrmsedge(18:27),'r.')
xlabel('\delta')
ylabel('rms edge strain')
legend('pre','post')

figure;
hold on
plot(exp4ddat(1:12),uxrmscent(1:12),'b.')
plot(exp4ddat(18:27),uxrmscent(18:27),'r.')
xlabel('\delta')
ylabel('rms cent strain')
legend('pre','post')



%%
%plotting as x = scan number
labels2 = cellstr(num2str([1:27]'));
figure;
dpts = 1:12;
hold on
plot(exp4ydat(dpts),'k*-','MarkerSize',10)
%plot((uxrmstot(dpts))*10000,'b*-','MarkerSize',10)
plot((uxrmscent(dpts))*10^4,'rx:','MarkerSize',10)
plot((uxrmsedge(dpts))*10^4,'go-','MarkerSize',10)
%plot(exp4ddat(dpts),'m*-','Markersize',10)
%plot(uxstdtot(dpts)*10000,'y*-','Markersize',10)
%title('Intermediate points ','FontSize',16)
xlabel('Scan Number','FontSize',16)
legend('OCV','rms cent x10^4','rms edge x10^4','rms edge (a.u.)','\delta in Li_{1-\delta}','std tot')
set(gca,'fontsize',16)
%text(1:length(dpts),exp4ydat(dpts),labels2(dpts), 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16)
hold off

%still organized by scan number, but look at rates of increase/decrease in
%center/edge
figure;
dpts = 2:12;
hold on
plot(exp4ydat(dpts-1),'k*-','MarkerSize',10)
plot((uxrmscent(dpts)-uxrmscent(dpts-1))*10000,'rx:','MarkerSize',10)
plot((uxrmsedge(dpts)-uxrmsedge(dpts-1))*10000,'g*-','MarkerSize',10)
legend('OCV','\Delta rms cent','\Delta rms edge')
set(gca,'fontsize',16)

%scatter plot for APS proposal
figure;
hold on
plot(exp4ddat,exp4ydat,'k*','MarkerSize',10)
xlabel('\delta','FontSize',16)
ylabel('OCV','FontSize',16)
set(gca,'fontsize',16)
title('Previous scan points','FontSize',16)


%including graphical representation of the extra cycles
labels = cellstr( num2str([1:17 zpts 18:27]') );
figure;
hold on
plot([exp4ydat(1:17) intpts exp4ydat(18:27)],'k*-','MarkerSize',10)
%plot([(uxrmstot(1:17))*10000 zpts (uxrmstot(18:27))*10000] ,'b*-','MarkerSize',10)
%plot([(uxrmscent(1:17))*10000 zpts (uxrmscent(18:27))*10000],'rx:','MarkerSize',10)
%plot([(uxrmsedge(1:17))*10000 zpts (uxrmsedge(18:27))*10000],'g*-','MarkerSize',10)
%plot([exp4ddat(1:17)' zpts exp4ddat(18:27)'],'y*-','Markersize',10)
%title('Experiment 4 first discharge/second charge','FontSize',16)
ylabel('OCV','FontSize',20)
%ylim([3 5]);
set(gca,'xtick',[])

%legend('OCV','rms tot','rms cent','rms edge','\delta in Li_{1-\delta}')
% text(1:36, [exp4ydat(1:17) intpts exp4ydat(18:27)]', labels, 'VerticalAlignment','bottom',...
%     'HorizontalAlignment','right','Fontsize',16)
text(1:17, exp4ydat(1:17), labels(1:17), 'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','Fontsize',20)
text(27:36,exp4ydat(18:27), labels(27:36), 'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','Fontsize',20)
set(gca,'fontsize',20)
hold off

%%
%Figure 3
%plotting as x = delta, still want to label points by their scan number
labels = cellstr(num2str([1:27]'));
dpts = 1:27;
figure;
hold on
plot(exp4ddat(dpts),uxrmstot(dpts)*10^4,'k*','MarkerSize',10)
%plot(exp4ddat(dpts),uxrmscent(dpts)*10000,'r*','MarkerSize',10)
%plot(exp4ddat(dpts),uxrmsedge(dpts)*10000,'b*','MarkerSize',10)
%plot(exp4ddat(18:27),uxrmstot(18:27),'bo','MarkerSize',10)
%plot(exp4ddat(13:17),uxrmstot(13:17),'gx','MarkerSize',10)
xlabel('\delta in Li_{1-\delta}','FontSize',20)
ylabel('RMS total strain x10^4','FontSize',20)
%ylabel('total RMS strain','FontSize',16)
%legend('Cent rms','Edge rms','Fontsize',16)
%title('Last discharge/second charge','FontSize',16)
set(gca,'fontsize',20)
%set(gca,'fontWeight','bold')
%text(exp4ddat(dpts), uxrmscent(dpts)*10^4, labels(dpts), 'VerticalAlignment','bottom',...
%    'HorizontalAlignment','right','Fontsize',16)
text(exp4ddat(dpts), uxrmstot(dpts)*10^4, labels(dpts), 'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','Fontsize',20,'color','blue')
hold off
%%
%plotting as y= OCV, x = delta, points labled by strain rms values

%looking at OCV as function of delta, plot labels as rms tot/std tot/max
%tot, etc
 labels2 = cellstr(num2str(round(uxmaxtot*1000000)'));
%labels2 = cellstr(num2str(round(uxstdtot*1000000)'));
%labels2(2) = cellstr('207');
labels1 = cellstr(num2str([1:27]'));
figure;
dpts = 1:27;
hold on
plot(exp4ddat(dpts),exp4ydat(dpts),'k*')
text(exp4ddat(dpts),exp4ydat(dpts),labels2(dpts), 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16,'color','red')
text(exp4ddat(dpts),exp4ydat(dpts),labels1(dpts), 'VerticalAlignment','top','HorizontalAlignment','left','Fontsize',16)
ylabel('OCV','fontsize',16)
xlabel('\delta in Li_{1-\delta}','FontSize',16)
title('Last charge - RMS tot','fontsize',16)
hold off
%%
%DQDV plot, Fig 2
%strain
labels2 = cellstr(num2str(round(uxrmstot*1000000)'));
%labels2 = cellstr(num2str(round(uxstdtot*1000000)'));
labels2(2) = cellstr('207');
uxrmstot(2) = 207/1000000;
%labels1 = cellstr(num2str([1:27]'));

%as function of voltage
% for i=1:12
%     %get corresponding dQ/dV y data point
%     if i<=6 %discharge one
%     ind = find(exp4ydat(i)>=dcdc1(:,1)); %first column is voltage
%     dqdv(i) = dcdc1(ind(1),2); %dQ/dV point corresponding to where I measured
%     elseif i>=6 && i<=12 %charge 2
%     ind = find(exp4ydat(i)<=dcc2(:,1));
%     dqdv(i) = dcc2(ind(1),2);
%     end   
% end

figure;
dpts = 1:12;
hold on
%text(exp4ydat(dpts),dcdc3(dqdv(:),2),labels2(dpts), 'VerticalAlignment','top','HorizontalAlignment','right','Fontsize',16,'color','red')
bar(exp4ydat(1:6),uxrmstot(1:6)*100,1)
bar(exp4ydat(7:12),uxrmstot(7:12)*100,1)
%plot(exp4ydat(dpts),dqdv,'r*') %voltage vs. dQ/dV CXDI points
plot(dcdc1(:,1),dcdc1(:,2),'k-','linewidth',6) %voltage vs. dQ/dV from before
plot(dcc2(:,1),dcc2(:,2),'r-','linewidth',6)
ylabel('dQ/dV','fontsize',16)
xlabel('OCV','fontsize',16)
title('Final Discharge - rmstot in red','fontsize',16)
hold off


%as function of delta



figure;
hold on;
% [ax,ch1,ch2]=plotyy(dvsVcharge(:,1),dvsVcharge(:,2),exp4ddat(7:12),uxrmstot(7:12),@plot,@bar); %voltage vs. dQ/dV from before
% axis(ax(1),[-.05 1 3 5.5]);
% axis(ax(2),[-.05 1 0 4e-4]);
bar(exp4ddat(2:6),uxrmstot(2:6)*10^4,1,'r')
bar(exp4ddat(7:12),uxrmstot(7:12)*10^4,1,'k')
 plot(dvsVcharge(:,1),dvsVcharge(:,2),'k-','linewidth',2) %voltage vs. dQ/dV from before
 plot(dvsVdischarge(:,1),dvsVdischarge(:,2),'r-','linewidth',2)
%[ax2,dc1,dc2]=plotyy(dvsVdischarge(:,1),dvsVdischarge(:,2),exp4ddat(2:6),uxrmstot(2:6),@plot,@bar); %voltage vs. dQ/dV from before
%axis(ax2(1),[-.05 1 3 5.5]);
%axis(ax2(2),[-.05 1 0 4e-4]);
ylabel('Voltage','fontsize',16)
ylim([0 5.5])
xlim([-0.05 1]);
xlabel('\delta in Li_{1-\delta}','FontSize',16)
legend('D.C. RMS total x10^5','C. RMS total x10^5')
title('First Curve','fontsize',16)
hold off




%%
%QUANT METRICS AS FUNCTION OF DELTA
deltaord = dscany;
deltax = dscanx;

for q=1:27

clear phasesp ampp ux;
%strain map
%pval = .03; %amplitude with connected pieces, play with this a bit
pval=0.2;
aeff2 = aeff(deltaord(q)); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

ampp = ampexp4(:,:,:,deltaord(q));
[xmax ymax zmax] = size(ampp);
phasesp = (phexp4(:,:,:,deltaord(q)).*(ampp>=pval))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan; %only want gradients where you have data
[ux uy uz] = gradient(phasesp,200,200,470); %use spacing of 20nm or 200 angstrom
clear uy uz;

%play around with this mask to try and get more reduction in surface than
%in center
%center
uxcent = (ux.*(ampp>=0.53));
uxcent(uxcent==0)=nan;
%edge gets 3/5 for surface dissolution 
% uxedge = (ux.*(ampp >= 0.19 & ampp<= 0.21));
% uxedge(uxedge==0)=nan;

%edge gets 4/5 for surface dissolution
uxedge = (ux.*(ampp >= 0.2 & ampp<= 0.25));
uxedge(uxedge==0)=nan;

%calculate the quant metrics
uxqtot(:,q) = reshape(ux,xmax*ymax*zmax,1);
uxqcent(:,q) = reshape(uxcent,xmax*ymax*zmax,1);
uxqedge(:,q) = reshape(uxedge,xmax*ymax*zmax,1);

uxnmean(q) = nanmean(uxqtot(:,q));
uxnmeancent(q) = nanmean(uxqcent(:,q));
uxnmeanedge(q) = nanmean(uxqedge(:,q));

%average of abs
uxavgtot(q) = nanmean(abs(uxqtot(:,q)));
uxavgcent(q) = nanmean(abs(uxqcent(:,q)));
uxavgedge(q) = nanmean(abs(uxqedge(:,q)));
%standard deviation
uxstdtot(q) = nanstd(uxqtot(:,q));
uxstdcent(q) = nanstd(uxqcent(:,q));
uxstdedge(q) = nanstd(uxqedge(:,q));
%maximums
uxmaxtot(q) = max(abs(uxqtot(:,q))); %this also might be useful
uxmaxcent(q) = max(abs(uxqcent(:,q)));
uxmaxedge(q) = max(abs(uxqedge(:,q)));
%rms
uxrmstot(q) = sqrt(nanmean(uxqtot(:,q).^2));
uxrmscent(q) = sqrt(nanmean(uxqcent(:,q).^2));
uxrmsedge(q) = sqrt(nanmean(uxqedge(:,q).^2));

end

labels = cellstr( num2str(dscany) );


%%

figure;
dpts = 1:27;
hold on
plot(deltax,(uxstdtot(dpts)*10000),'b*-','MarkerSize',10)
%plot(deltax,(uxrmscent(dpts)),'rx:','MarkerSize',10)
%plot(deltax,(uxrmsedge(dpts)),'g*-','MarkerSize',10)
title('Experiment 4 all data','FontSize',16)
xlabel('\delta in Li_{1-\delta}','FontSize',16)
legend('rmstot','rms cent','rms edge')
text(deltax, uxstdtot*10000, labels, 'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','Fontsize',16,'color','red')
set(gca,'fontsize',16)
hold off

%%
labels2 = cellstr(num2str([1:27]'));
figure;
hold on
plot(1:27,exp4ydat,'k*-')
plot(1:27,exp4ddat,'b*-')
legend('OCV','\delta')
text(1:27,exp4ydat,labels2, 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16)
ylabel('OCV','fontsize',16)
set(gca,'fontsize',16)
hold off


%CAPACITY FADE CORE
%going to need a different one for the higher z resolution guys
avgamp2 = mean(ampexp4(:,:,:,[15:17 22 27]),4);
%and also don't forget to change the cmask stuff
capdat = [1 6 12 13 14 15 16 17 18 22 27];

%1,6,12-14 use avgamp, 15 16 17 22 27 use avgamp2, 18 uses its own

for q=1:11
clear phasesp ampp ux;
%strain map
%pval = .03; %amplitude with connected pieces, play with this a bit
pval=0.1;
aeff2 = aeff(capdat(q)); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

ampp = ampexp4(:,:,:,capdat(q));
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
phasesp = (phexp4(:,:,:,capdat(q)).*(ampp>=pval))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,200,200,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
ux(isnan(ux)==1) = 0;

figure;
hold on
cdata = smooth3(ux,'box',[3 3 1]);
%decide if you wanna use 0.8 or 0.5 here
p = patch(isosurface(x,y,z,ampp,0.5));
isonormals(x,y,z,ampp,p);
isocolors(x,y,z,cdata,p);
set(p,'FaceColor','interp','EdgeColor','none')
daspect([1 1 1]);axis tight
camlight(230,250); 
lighting GOURAUD; 
%title('4.85 2nd Charge ','FontSize',16)
%colorbar;
set(gca, 'ydir','reverse')
view(3); 
axis off;
caxis([-5e-4 5e-4]); 
set(gca,'fontsize',16)
hold off

end





%testing out idea of finding where the region of highest strain is
for q=1:12
pval2 = 0.2; %doesn't matter so much, use .03 or 0.1
aeff2 = aeff(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

if szvec(q) == 64
    ampp = avgamp;
elseif szvec(q) == 96
    ampp = avgamp2;
elseif szvec(q) == 80
    ampp = avgamp3;
end

phasesp = (phexp4(:,:,:,q).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
[ux uy uz] = gradient(phasesp,200,200,470);
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
cdata = smooth3(ux,'box',[3 3 1]);

[r,c,v] = ind2sub(size(cdata),find(abs(cdata) >= 0.9* max(max(max(abs(cdata))))));

figure;
p = patch(isosurface(x,y,z,ampp,0.2));
hold on
set(p,'faceAlpha',0.2)
plot3(c,r,v,'r*','Markersize',30)
view(3); 
axis off;
%text(c,r,v, labels, 'VerticalAlignment','bottom',...
%    'HorizontalAlignment','right','Fontsize',30,...
%    'BackgroundColor',[.7 .9 .7],'EdgeColor','red','LineWidth',3);

clear r c v;
end


%%
labels = cellstr(num2str([1:6]'));
figure;
p = patch(isosurface(x,y,z,ampp,0.2));
hold on
set(p,'faceAlpha',0.3)
plot3(c,r,v,'r*','Markersize',30)
view(3); 
axis off;
text(c,r,v, labels, 'VerticalAlignment','bottom',...
    'HorizontalAlignment','right','Fontsize',30,...
    'BackgroundColor',[.7 .9 .7],'EdgeColor','red','LineWidth',3);



%DELTA ORGANIZED GRAPHS
deltaord = dscany;
deltax = dscanx;

%want to go in order of decreasing delta/increasing Li
%can try changing pval, ampp, where you slice
for q=1:27

pval2 = 0.1; %doesn't matter so much, use .03 or 0.1
aeff2 = aeff(deltaord(q)); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

if szvec(deltaord(q)) == 64
    ampp = avgamp;
elseif szvec(deltaord(q)) == 96
    ampp = avgamp2;
elseif szvec(deltaord(q)) == 80
    ampp = avgamp3;
end

phasesp = (phexp4(:,:,:,deltaord(q)).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,200,200,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
uxnans = ux;
cdata = smooth3(ux,'box',[3 3 1]);
ux(isnan(ux)==1) = 0;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);


figure;
imagesc(cdata(:,:,48));
%caxis([-1e-4 1e-4]) %used to use -5e-5 5e-4
zoom(4)
axis off;

%bin data
uxnans(uxnans > 1e-4) = 1e-4; 
uxnans(uxnans < -1e-4) = 1e-4;
ux(ux > 1e-4) = 1e-4; 
ux(ux < -1e-4) = 1e-4;

% figure;
% imagesc(ux(:,:,48));
% %caxis([-5e-4 5e-4]) %used to use -5e-5 5e-4
% zoom(4)
% axis off;

uxnansedit = reshape(uxnans(:,:,48),xlim*ylim,1);

%it looks like this way will work
varplot(q) = nanstd(uxnansedit);
end

labels2 = cellstr(num2str(deltaord));
figure;
hold on
plot(deltax,varplot,'k*-')
text(deltax,varplot,labels2, 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16)
xlabel('delta in Li_{1-\delta}','fontsize',16)
ylabel('Variance after binning - high implies small stripes','fontsize',16)
hold off

labels2 = cellstr(num2str(deltaord));
figure;
hold on
plot(licontent,stripewidth,'k*-')
text(licontent,stripewidth,labels2, 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16)
xlabel('Li content','fontsize',16)
ylabel('1/variance - measure of stripe width','fontsize',16)
hold off

%histogram investigation, two phase is 10/22
q=4;
figure;
hist(uxqtot(:,q),[-6e-4:.5e-4:6e-4])
xlabel('values of strain','fontsize',16)
ylabel('number of occurances','fontsize',16)
%title('two phase histogram - first cycle','fontsize',16)


labels2 = cellstr(num2str([1:27]'));
figure;
hold on
plot(1:27,exp4ddat,'k*-')
text(1:27,exp4ddat,labels2, 'VerticalAlignment','bottom','HorizontalAlignment','right','Fontsize',16)
ylabel('delta in Li_{1-\delta}','fontsize',16)
hold off


%plotting variance as a function of radial position

for q=1:27

clear phasesp ampp ux;
%strain map, was using 0.1 before
pval=0.2;
aeff2 = aeff(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

ampp = ampexp4(:,:,:,q);
[xmax ymax zmax] = size(ampp);
phasesp = (phexp4(:,:,:,q).*(ampp>=pval))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan; %only want gradients where you have data
[ux uy uz] = gradient(phasesp,200,200,470); %use spacing of 20nm or 200 angstrom
clear uy uz;

%break up into 6 regions
ux1 = reshape((ux.*(ampp >= 0.2 & ampp<= 0.25)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(1,q) = sqrt(nanmean(ux1.^2));
clear ux1;

ux1 = reshape((ux.*(ampp >= 0.26 & ampp<= 0.3)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(2,q) = sqrt(nanmean(ux1.^2));
clear ux1;

ux1 = reshape((ux.*(ampp >= 0.3 & ampp<= 0.35)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(3,q) = sqrt(nanmean(ux1.^2));
clear ux1;

ux1 = reshape((ux.*(ampp >= 0.35 & ampp<= 0.45)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(4,q) = sqrt(nanmean(ux1.^2));
clear ux1;

ux1 = reshape((ux.*(ampp >= 0.45 & ampp<= 0.65)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(5,q) = sqrt(nanmean(ux1.^2));
clear ux1;


ux1 = reshape((ux.*(ampp >= 0.65 & ampp<= 1)),xmax*ymax*zmax,1);
ux1(ux1==0)=nan;
ux1 = snip(ux1,nan);
nnz(ux1)
ux1rms(6,q) = sqrt(nanmean(ux1.^2));
clear ux1;



end


figure;
hold on
plot(ux1rms)
xlabel('distance from edge','fontsize',16)
ylabel('RMS strain','fontsize',16)
hold off


%showing distribution of strains around the effective lattice constant
%pretty cool, wait for feedback from oleg/shirley but I like it

x=1:12;
y=aeff4(x);

%this will work beautifully!!
hold on
for i=2:12
h=plot(x(i),y(i),'ro');
set(h,{'markers'},{uxrmscent(i)*200000}) 

h=plot(x(i),y(i),'bo');
set(h,{'markers'},{uxrmsedge(i)*200000}) 

plot(x(i),y(i),'k.','markersize',10);
end


%make a rotating movie of the particle
% Prepare the new file.
    vidObj = VideoWriter('part_rot.avi');
    open(vidObj);
 
% Create an animation.
view(3);
pos = get(gcf,'position');
 
 for i=0:10:360
        view(0+i,30)
        set(gcf,'position',pos)
       % Write each frame to the file.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
 end
  
    % Close the file.
    close(vidObj);
		
%% Make my own PXRD graph

hold on
%charge data
plot(exp4ddat(7:12),aeff4(7:12),'bo') %my individual particle data
plot(exp4ddat(1:6),aeff4(1:6),'k+') %my individual particle data
plot(josh_delta,josh_lat1,'rx') %josh phase 1
plot(josh_delta(end-3:end),josh_lat2,'rx') %josh phase 1 %josh phase 2
plot(josh_delta(end-2:end),josh_lat3,'rx') %josh phase 1 %josh phaes 3
%discharge data
hold on
plot(josh_delta_d,josh_lat1_d,'kx') %josh phase 1
plot(josh_delta_d(1:2),josh_lat2_d,'kx') %josh phase 1 %josh phase 2
plot(josh_delta_d(1),josh_lat3_d,'kx') %josh phase 1 %josh phaes 3
plot(exp4ddat(1:6),aeff4(1:6),'r+') %my individual particle data
xlabel('\delta in Li_{1-\delta}')
ylabel('Lattice Constant (Angstrom)')






%%
%do the strain energy calculation
Y = 150e9; %200 GPa Youngs modulus
K = 130e9; %120 GPa Bulk modulus

% C11=220e9; %used for poisson ratio
% syms x;
% solve(C11==(1-x)*Y/((1+x)*(1-2*x)),x)

n = 0.327; %poisson ratio I found using solve

G = Y/(2*(1+n));
I = (2*n*G)/(1-2*n);

%grams of Li in particle

gli = 8*(sqrt(2)/3)*((400e-9)^3)/((8.1e-10)^3)*6.94/(6.022e23);



capty = [72 23 22 20 8 0 10 21 75 90 110 140]; %mAh/g

kap = 5.02e-10; %J/m

for q=[7:12]
pval2 = 0.1; %doesn't matter so much, use .03 or 0.1
aeff2 = aeff4(q); %this changes for every scan, put in Angstrom
phsf =  2*pi*sqrt(3)/aeff2; %in Angstroms now

if szvec(q) == 64
    ampp = avgamp;
elseif szvec(q) == 96
    ampp = avgamp2;
elseif szvec(q) == 80
    ampp = avgamp3;
end

phasesp = (phexp4(:,:,:,q).*(ampp>=pval2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
uxnans = ux;

[uxx uyy uzz] = gradient(ux,16e-9,16e-9,47e-9);

%ux(isnan(ux)==0) = (ux(isnan(ux)==0) + avgstrainabs(q)); %put in the average abs strain

%sum of the squared values times the voxel volume
%Estrain(q) = (2*G+3*I)/2*nansum(nansum(nansum(ux.^2)))*((16e-9)*(16e-9)*(47e-9));

%do it adding in the minimum strain to everything for q=1:6
%need to take region and average to get the minimum
cdata = smooth3(ux,'box',[5 5 3]);
%ux(isnan(ux)==0) = ux(isnan(ux)==0) + abs(min(min(min(cdata))));
%do it subtracting the max going the other way q=7:12
ux(isnan(ux)==0) = ux(isnan(ux)==0) - max(max(max(cdata)));

Estrain(q) = (2*G+3*I)/2*nansum(nansum(nansum(ux.^2)))*((16e-9)*(16e-9)*(47e-9));


%Evolt(q) = (1-exp4ddat(q))*exp4ydat(q)*capty(q)*(1/1000)*3600*gli;

Ekap(q) = kap*nansum(nansum(nansum(uxx.^2+uyy.^2+uzz.^2)))*((16e-9)*(16e-9)*(47e-9));

end

fdc=1:6;
figure;
plot(exp4ddat(fdc),Estrain(fdc),'kd') %black diamonds for discharge
%xlabel('\delta')
%ylabel('Mechanical Energy')
%title('First Charge')
xlim([0 1]);
%ylim([0 1.5e-12]);
saveas(1,'J4FDC.eps','psc2')

sc=7:12;
figure;
plot(exp4ddat(sc),Estrain(sc),'bs') %blue square for charge
%xlabel('\delta')
%ylabel('Mechanical Energy')
%title('First Charge')
xlim([0 1]);
%ylim([0 1.5e-12]);
saveas(1,'J4SC.eps','psc2')

figure; hold on;
plot(exp4ddat(sc),Estrain(sc),'bs') %blue square for charge
plot(exp4ddat(fdc),Estrain(fdc),'kd') %black diamonds for discharge
%xlabel('\delta')
%ylabel('Mechanical Energy')
%title('First Charge')
xlim([0 1]);
ylim([0 2e-14]);
saveas(1,'both.eps','psc2')




%% Playing with multiple contour slice visualization

load('/Users/andrewulvestad/Documents/exp4Oct25.mat')
avgamp = mean(ampexp4(:,:,:,1:14),4); %use first 12 bc all have 64 z dim
ampp = avgamp;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
phsf =  2*pi*sqrt(3)/8.14;
set(gcf, 'InvertHardCopy', 'off');

for q=[12]

phasesp = (phexp4(:,:,:,q).*(ampp>=0.2))/phsf; %this is now displacement in angstrom
phasesp(phasesp==0)=nan;
[ux uy uz] = gradient(phasesp,160,160,470); %use spacing of 20nm or 200 angstrom
clear uy uz;
ux(isnan(ux)==1) = 0;
[ylim xlim zlim] = size(ampp);
[x y z] = meshgrid(1:xlim,1:ylim,1:zlim);
cdata = smooth3(ux,'box',[3 3 1]);


figure;
%this is a really great viz tool
h = contourslice(x, y, z, cdata, [], [], [96/2-2 96/2-1 96/2 96/2+1 96/2+2],256);
camva(24); camproj perspective;
campos([-3, -15, 5]);
set(gcf, 'Color', [.3, .3, .3], 'renderer', 'zbuffer');
set(gca, 'Color', 'black', 'XColor', 'black', ...
'YColor', 'black', 'ZColor', 'black');
box off;
view(-38,16)

caxis([-4e-4 4e-4]);
zoom(3.5);
axis off;
%set(gca,'color','k')
set(h, 'LineWidth', 6);
set(gcf,'color','w')
set(gcf, 'InvertHardCopy', 'off');
print -deps2c output.eps
end








