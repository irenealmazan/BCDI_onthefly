%% on the fly for loading tiff stack, remove aliens, put into data and run it quickly
clear all; close all; clc; %initialize_plot_PRL;

% load in tiff stack

cd ~/Desktop/
curdir = pwd; %start on desktop

scannum2=[874 932];
viewdat =0;
cnt=1;
wantin = 0;
remalienfromprev = 0;
%
for   i=1:length(scannum2)
    close all;
    scannum = scannum2(i);
    
    %b
    datdir = sprintf('/Volumes/34idc-acq/2016/You1116/You1116c_S%04d/',scannum);
		
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
    
if nfilesmax>25 
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

%alien removal step      
% for gw=3:2:nfilesmax
% 			imagesc(log(data(:,:,gw)))
%             title(['zsl ' num2str(gw) ' diff data for scannum' num2str(scannum)])
%     %        imagesc(data(:,:,g))
% 			%pause(0.2)
%             
%             result = input('require alien removal?');
%  
% % 
%                 if result>=1
%                     for cvg = 1:result
%                 disp('click on center of alien, surrounding 20 pixels will be zerod');
%                 [xc yc] = ginput(1);
%                 xc = round(xc);
%                 yc = round(yc);
%                 data(abs(yc-10:yc+10)+1,abs(xc-10:xc+10)+1,abs(gw-2:gw+2))=0;
%                 
%                 xcstore(cnt) = xc;
%                 ycstore(cnt) = yc;
%                 zcstore(cnt) = gw;
%                 
%                 cnt = cnt+1;
%                 
%                 clear xc yc;
%                     end
%                 end
% 
% end

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
     
     % compute the lattice constant
medfiltuse=0;
nofilt=0;
smth3=1;
    
cd('/Volumes/34idc-acq/2016/You1116/')

[specscan, errors] = openspec('You1116c.spec',scannum);

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
th_ofset = 0.1;
%th_ofset = 0;-.01%62;

a_eff(scannum) = sqrt(3)*1.37761/(2*sind((eff_2theta(scannum)+th_ofset)/2));

widthx(scannum) = pkwidth_x;
widthy(scannum) = pkwidth_y;

totali(scannum) = sum(data(:));
 
 
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
a_eff(a_eff~=0)
figure; plot(a_eff(a_eff~=0));
