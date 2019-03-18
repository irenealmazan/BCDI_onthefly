%% on the fly for loading tiff stack, remove aliens, put into data and run it quickly
clear all; close all; clc; %initialize_plot_PRL;

% load in tiff stack

cd ~/Desktop/
curdir = pwd; %start on desktop

scannum2=[123];
viewdat =0;
cnt=1;
wantin = 0;
remalienfromprev = 0;
%
for   i=1:length(scannum2)
    close all;
    scannum = scannum2(i);
    
    datdir = sprintf('/Volumes/34idc-acq/2016/Hruszkewycz1216/Hruszkewycz1216a_S%04d/',scannum);
		
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
 %	 save(sname,'data');
     
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
play_data(data);
