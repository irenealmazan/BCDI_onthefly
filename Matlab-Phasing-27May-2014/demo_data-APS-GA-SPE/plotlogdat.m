%get data and save it
clear all; close all; clc;

cd ~/Desktop/

plotlogdat = 1;
centdat=1;
%scannum2=[2857, 2881	2898	2926, 2930	2958, 2962	2982	2998	3014	3030	3046	3070	3086	3134	3146	3158	3170:3195];

scannum2=[1460 1487 1495];

%load it up 
figure;

 initialize_plot_PRL;
for i=1:length(scannum2)
    
    scannum = scannum2(i);
    
    try

    
 dname = sprintf('S%04d.mat',scannum);
 
 load(['~/Documents/ulvestad416b_dsets/' dname])

 
  if centdat == 1
 data=center_array(data);
 disp('Centering data....')
 end

  imagesc(log10(data(:,:,ceil(end/2)))); axis equal; colorbar;
  title(['log diff data for ' num2str(scannum)])
  cd ~/Desktop/
    if plotlogdat == 1
   
  imagesc(log10(data(:,:,ceil(end/2)))); axis equal; axis off; 
    saveas(gcf,[ num2str(scannum) '_nolabels.pdf'],'pdf')
    end
 
    catch me
    end
end