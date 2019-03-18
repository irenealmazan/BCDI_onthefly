%% look at the data quality judged as the integrated number of photons
cd ~/Desktop/
%make this standalone can just grab the centers of the data
clear all; close all; clc;
scannum2=[283:291];

cnt = 1;
mindata=2;
centdat=1;
pldat=0;
initpad=1;
plotcent =1;
 %initialize_plot_PRL;
for i=scannum2
    
 dname = sprintf('S%04d.mat',i);
 
 try 
 
 load(['~/Desktop/' dname])
 
 ind=( data < mindata );
 data(ind)=0;
 
%  if initpad ==1
%  data = init_pad(data,[10 10 10 10 0 0]);
%  data=fft_padjesse(data,[[1 1],1]);
%  end
 
 
%  if centdat == 1
%  data=center_array(data);
%  disp('Centering data....')
%  end
%  
%  if pldat == 1
%      play_data(data)
%  end
%  
%  data = sqrt(data);
 
  [p1 p2 p3] = ind2sub(size(data),find(data>=1000));
   %condition of max being greater than 150% of n.n.
   %length(p1)
   for qq=1:length(p1)
   %disp('checking for saturated pixels')
   nn1 = data(p1(qq)+1,p2(qq),p3(qq));
   nn2 = data(p1(qq)-1,p2(qq),p3(qq));
   nn3 = data(p1(qq),p2(qq)+1,p3(qq));
   nn4 = data(p1(qq),p2(qq)-1,p3(qq));
   mval = data(p1(qq),p2(qq),p3(qq));
   
   if mval>500 && mval>=1.5*nn1 && mval>=1.5*nn2 && mval>=1.5*nn3 && mval>=1.5*nn4
       disp('saturated pixel found with value')
       mval
       disp('changing value to')
       newval = (1/4)*(nn1+nn2+nn3+nn4);
 %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
 %      pause
       data(p1(qq),p2(qq),p3(qq))=newval;
 %      scannum
 %      close(gcf)
   end
   end

cd('/Volumes/34idc-acq/2016/You1116/')

[specscan, errors] = openspec('You1116b.spec',i);

cd('~/Desktop/')

 
%  if plotcent ==1
%   for zz=1:size(data,3)
%       dfilt(:,:,zz) = medfilt2(data(:,:,zz));
%   end
%      figure; imagesc(log10(data(:,:,end/2-5))); axis equal; colorbar;  drawnow;
%  end
 
 
 totali(:,cnt) = squeeze(sum(squeeze(sum(data,1)),1));
 
 totalmon(:,cnt) =  specscan.data(21,:);
 
 totalinorm(:,cnt) = totali(:,cnt)./totalmon(:,cnt);
 
 Estore(cnt) = specscan.motor_positions(30);
 
 %data = data/max(data(:));
 
 sumdat2d(:,:,cnt) = sum(data,3);

 data(data==0)=nan;
 dstore(:,cnt) = data(:);

clear scannum;


clear data;
snumstore(cnt) = i;
cnt = cnt+1;
 catch me
     me.message
 end

end

figure; subplot(2,2,1); plot(totalinorm(:,1),'r'); subplot(2,2,2); plot(totalinorm(:,2),'k')
%%
% correlation analysis
[r,p] = corrcoef(dstore,'rows','pairwise');

figure; imagesc(r); drawnow;

[rlog,p] = corrcoef(log10(dstore),'rows','pairwise');

figure; imagesc(rlog); drawnow;

figure; plot(totali,'b-.')

%% subplot the data
ndset = length(snumstore);
for i=1:ndset
    subplot(ceil(sqrt(ndset)),ceil(sqrt(ndset)),i); imagesc(log10(sumdat2d(:,:,i)))
end



