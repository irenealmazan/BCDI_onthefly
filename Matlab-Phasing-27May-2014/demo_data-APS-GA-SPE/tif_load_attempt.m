%stuff for opening all the tif files from experiment


%Staff1212_S120_00000

% scannum = 120;
% 
% datdir = ['/Users/a4894z/Desktop/S' num2str(scannum) '/'];
% 
% nfiles = 51;
% 
% sl = imread([datdir 'Staff1212_S' num2str(scannum) '_000']
% 
% 
% params.files={['Staff1212_S' num2str(scannum) '_00000.tif'],['Staff1212_S' num2str(scannum) '_00050.tif']} ;
% params.back={['Staff1212_S' num2str(scannum) '_bg00051.tif']};

%/Users/a4894z/Desktop/

%first load the tifs and put into a matlab stack
%clear all; clc; close all;

i=0;
clear files2;
curdir = pwd; %start on desktop

% [702:715 717:731 733:987] november particle 4

for scannum = [750:754 840:850]
    datdir = sprintf('/Users/a4894z/Documents/oleg1113/oleg1113_S%04d/',scannum);
    datdir2 = sprintf('/Users/a4894z/Documents/oleg1113/oleg1113_S%04d',scannum);
    cd(datdir);
    clear listOftifs;
		clear dir;
    listOftifs = dir('*.tif');
    nfilesmax = numel(listOftifs)
    data = zeros(256,256,nfilesmax);
    cd(curdir);
    
    for nfiles = 0:nfilesmax-1
        filename = sprintf('oleg1113_S%04d_%05d.tif',scannum,nfiles);
         I = imread([datdir filename]);
         data(:,:,nfiles+1) = I;
        %fprintf('Loading file %s \n',[datdir filename])
    end
    sname = sprintf('S%04d',scannum);
    i=i+1;
    
    
		
		
		%get rid of the saturated pixels somehow
  % [p1 p2 p3] = ind2sub(size(data),find(data>=0.9*max(max(max(data)))));
  [p1 p2 p3] = ind2sub(size(data),find(data>=1000));
   %condition of max being greater than 150% of n.n.
   %length(p1)
   for qq=1:length(p1)
   
   nn1 = data(p1(qq)+1,p2(qq),p3(qq));
   nn2 = data(p1(qq)-1,p2(qq),p3(qq));
   nn3 = data(p1(qq),p2(qq)+1,p3(qq));
   nn4 = data(p1(qq),p2(qq)-1,p3(qq));
   mval = data(p1(qq),p2(qq),p3(qq));
   
   if mval>500 && mval>=1.5*nn1 && mval>=1.5*nn2 && mval>=1.5*nn3 && mval>=1.5*nn4
       disp('saturated pixel found with value')
       mval
       disp('changing value to')
       newval = (1/4)*(nn1+nn2+nn3+nn4)
 %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
 %      pause
       data(p1(qq),p2(qq),p3(qq))=newval;
 %      scannum
 %      close(gcf)
   end
   end
   
    %check the slices after removal
    
     %look at the slices individually
% 		figure;
% 	for g=1:nfilesmax
% 			imagesc(log(data(:,:,g)))
%     %        imagesc(data(:,:,g))
% 			pause(0.4)
%             %pause;
%   end
%     close gcf;
    
    
    
    %finally save the data
    save(sname,'data');
  %  sname = [sname '.mat'];
 %   files2{i}= sname;
		
   
    
%       % Prepare the new file.
%       vidObj = SpeWriter('test.spe');
%       open(vidObj);
%     
%       % Create an animation.
%     
%       for k = 1:nfilesmax
%          imshow(data(:,:,k)) 
%          % Write each frame to the file.
%          currFrame = getframe;
%          writeVideo(vidObj,currFrame);
%       end
%     
%       % Close the file.
%       close(vidObj);
    
% find the coordinates of the maximum        
%[p11 p21 p31] = ind2sub(size(data),find(data==max(max(max(data)))));

%store the coordinates
%bkpkcents(scannum,1:2) = [p11(1) p21(1)];



%look at the summed image
% imagesc(log(sum(data,3))); pause(0.5)
%  
%fit a gaussian to the data
% ypts = sum(sum(data,3),1);
% xpts = 1:length(ypts);
% f=fit(xpts.',ypts.','gauss1');
% 
% %find the location of the maximum
% pkloc(scannum) = find(ypts==max(ypts));
% 
% %parameters from the gaussian
% afits(scannum) = f.a1;
% bfits(scannum) = f.b1;
% cfits(scannum) = f.c1;
% 
% lnespl(:,scannum)=ypts;


clear data;
end

% hold on
% for zz=1152:1226
%     %plot3(1:length(ypts),lnespl(:,zz),ones(1,length(ypts))*zz)
% 		plot(1:length(ypts),lnespl(:,zz))
% 		pause(0.25)
% end


%%
%for normal loading of scans

data = zeros(256,256,31);

for scannum = 411:420
    datdir = sprintf('/Users/andrewulvestad/Documents/oleg1113/oleg1113_S%04d/',scannum);
    %datdir = sprintf('/Users/a4894z/Documents/oleg1113/oleg1113_S%04d/',scannum);
    for nfiles = 0:30
        filename = sprintf('oleg1113_S%04d_%05d.tif',scannum,nfiles);
         I = imread([datdir filename]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir filename])
    end
    sname = sprintf('S%04d',scannum);
    save(sname,'data');
    data = zeros(256,256,31);
end






%for the DS ring analysis





data = zeros(1340,300,566); %rows, cols

for scannum = 11:11
    datdir = sprintf('oleg1113b_S%04d/',scannum);
    for nfiles = 0:565
        filename = sprintf('oleg1113b_S%04d_%05d.tif',scannum,nfiles);
         I = imread([datdir filename]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir filename])
    end
end

bkg = zeros(1340,300,10);

for scannum = 11:11
    datdir = sprintf('oleg1113b_S%04d/',scannum);
    for nfiles = 569:577
        filename = sprintf('oleg1113b_S%04d_%05d.tif',scannum,nfiles);
         I = imread([datdir filename]);
         bkg(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir filename])
    end
end


bkg(:,:,1:568)=[];
