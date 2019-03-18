%% look at the data quality judged as the integrated number of photons

%make this standalone can just grab the centers of the data
clear all; close all; clc;
%p1
scannum2=[2013];

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
 
 load(['~/Documents/Hoydoo216bdsets/' dname])
 
 ind=( data < mindata );
 data(ind)=0;
 
 if initpad ==1
 data = init_pad(data,[10 10 10 10 0 0]);
 data=fft_padjesse(data,[[1 1],1]);
 end
 
 
 if centdat == 1
 data=center_array(data);
 disp('Centering data....')
 end
 
 if pldat == 1
     play_data(data)
 end
 
 data = sqrt(data);
 
 
 if plotcent ==1
  for zz=1:size(data,3)
      dfilt(:,:,zz) = medfilt2(data(:,:,zz));
  end
     figure; imagesc(log10(data(:,:,end/2-5))); axis equal; colorbar;  drawnow;
 end
 
 
 totali(cnt) = sum(data(:));
 
 data = data/max(data(:));
 
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



