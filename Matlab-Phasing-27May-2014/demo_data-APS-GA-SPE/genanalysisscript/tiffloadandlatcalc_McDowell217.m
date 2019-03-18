%% on the fly for loading tiff stack, remove aliens, put into data and run it quickly
clear all; close all; clc; %initialize_plot_PRL;

% load in tiff stack

cd ~/Desktop/
curdir = pwd; %start on desktop

scannum2=[190:260];
viewdat =0;
cnt=1;
wantin = 0;
remalienfromprev = 0;
%
for   i=1:length(scannum2)
    close all;
    scannum = scannum2(i);
    
    %b
    %datdir = sprintf('/net/s34data/export/34idc-acq/2017/McDowell217/McDowell217a_S%04d/',scannum);
     datdir = sprintf('~/Documents/McDowell217/McDowell217a_S%04d/',scannum);
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
    
if nfilesmax<100 
    data = zeros(256,256,nfilesmax);
    cd(curdir);
    
    for nfiles = 0:nfilesmax-1
         I = imread([datdir listOftifs(nfiles+1).name]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir listOftifs(nfiles+1).name])
    end
    sname = sprintf('S%04d',scannum);
       
     
     % compute the lattice constant
medfiltuse=0;
nofilt=0;
smth3=1;
    
cd('~/Documents/McDowell217/')

[specscan, errors] = openspec('McDowell217a.spec',scannum);

cd ~/Desktop/

 Delta = specscan.motor_positions(1);
 Gamma = specscan.motor_positions(6);
 camdist = specscan.motor_positions(27)/1000 % in meters

  if medfiltuse==1
  %median filter first
  for zz=1:size(data,3)
      dfilt(:,:,zz) = medfilt2(data(:,:,zz));
  end
  
  %find slice with max
  [p1 p2 p3] = ind2sub(size(dfilt),find(dfilt==max(max(max(dfilt)))));
  elseif nofilt==1
      [p1 p2 p3] = ind2sub(size(data),find(data==max(max(max(data)))));
  elseif smth3==1
      dfilt = smooth3(data,'box',[3 3 1]);
      [p1 p2 p3] = ind2sub(size(data),find(dfilt==max(max(max(dfilt)))));
      disp('used smooth3 data')
  end
  
  ypts = sum(data(:,:,p3(1)),1);

xpts = 1:length(ypts);

% 
 %find the location of the maximum
% pkloc_x = find(ypts==max(ypts));

%use a gauss fit instead
f=fit(xpts.',ypts.','gauss1');
pkloc_x = f.b1;
pkwidth_x = f.c1;


% ypts2 = sum(sum(data,3),2);
ypts2 = sum(data(:,:,p3(1)),2);

 %pkloc_y = find(ypts2==max(ypts2));
 
 
 fy=fit(xpts.',ypts2,'gauss1');
 
 pkloc_y = fy.b1;
 pkwidth_y = fy.c1;
 
 
% 
 bkpkcents(scannum,1:2) = [pkloc_y(1) pkloc_x(1)];
 
 
xposCCD = bkpkcents(scannum,2); %colums
yposCCD = bkpkcents(scannum,1); %rows




%pix_delta = 128-xposCCD;
pix_delta = xposCCD-128;

pix_gamma = 128-yposCCD;
%pix_gamma = -128+yposCCD;

%cam dist and pixel size
add_delta = atand(((55e-6)*pix_delta)/(camdist));
add_gamma = atand(((55e-6)*pix_gamma)/(camdist));

eff_delta = add_delta+Delta; %these are in degrees

%add_gamma = 0;

eff_gamma = add_gamma+Gamma;

eff_2theta(scannum) = acosd(sind(90-eff_gamma)*cosd(eff_delta)); %result in degrees

%th_ofset = .0707;
th_ofset = -.5%.1;
%th_ofset = 0;-.01%62;

a_eff1(scannum) = sqrt(3)*1.37761/(2*sind((eff_2theta(scannum)+th_ofset)/2));

widthx(scannum) = pkwidth_x;
widthy(scannum) = pkwidth_y;

totali(scannum) = sum(data(:));
 
 
elseif nfilesmax>100
    cd(curdir);
    
    for nfiles = 0:nfilesmax-1
         I = imread([datdir listOftifs(nfiles+1).name]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir listOftifs(nfiles+1).name])
    end
    sname = sprintf('S%04d',scannum);
       
     
     % compute the lattice constant
medfiltuse=0;
nofilt=0;
smth3=1;
    
cd('~/Documents/McDowell217/')

[specscan, errors] = openspec('McDowell217a.spec',scannum);

cd ~/Desktop/

 Delta = specscan.motor_positions(1);
 Gamma = specscan.motor_positions(6);
 camdist = specscan.motor_positions(27)/1000 % in meters
    
for lne=1:nfilesmax    
ypts = sum(data(:,:,lne),1);

xpts = 1:length(ypts);

% 
 %find the location of the maximum
 pkloc_x = find(ypts==max(ypts));

%use a gauss fit instead
% f=fit(xpts.',ypts.','gauss1');
% pkloc_x = f.b1;
% pkwidth_x = f.c1;


% ypts2 = sum(sum(data,3),2);
ypts2 = sum(data(:,:,lne),2);

pkloc_y = find(ypts2==max(ypts2));
 
 
%  fy=fit(xpts.',ypts2,'gauss1');
%  
%  pkloc_y = fy.b1;
%  pkwidth_y = fy.c1;
 
 
% 
 bkpkcents(lne,1:2) = [pkloc_y(1) pkloc_x(1)];
 
 
xposCCD = bkpkcents(lne,2); %colums
yposCCD = bkpkcents(lne,1); %rows




%pix_delta = 128-xposCCD;
pix_delta = xposCCD-128;

pix_gamma = 128-yposCCD;
%pix_gamma = -128+yposCCD;

%cam dist and pixel size
add_delta = atand(((55e-6)*pix_delta)/(camdist));
add_gamma = atand(((55e-6)*pix_gamma)/(camdist));

eff_delta = add_delta+Delta; %these are in degrees

eff_gamma = add_gamma+Gamma;

eff_2theta(lne) = acosd(sind(90-eff_gamma)*cosd(eff_delta)); %result in degrees

%th_ofset = .0707;
th_ofset = -0.2;
%th_ofset = 0;-.01%62;

if max(ypts2)>60
a_eff2(lne) = sqrt(3)*1.37761/(2*sind((eff_2theta(lne)+th_ofset)/2));
end

end
    
    
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
a_eff1(a_eff1~=0)
figure; plot(a_eff1(a_eff1~=0),'k.-','MarkerSize',20);

%% convert to concentration of Ni

xNi = (3.89-a_eff1(a_eff1~=0))./.37;