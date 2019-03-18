function [data] = loadtiffdirectnewdet(path,scannum,params)

%load tif directly from folder
viewdat =0;
cnt=1;
wantin = 0;
remalienfromprev = 0;
padme=1;
centdat=1;
mindat=params.min_data;
curdir = params.curdir;

for   i=1
    close all;

    
    datdir = sprintf(['' path 'S%04d/'],scannum);
		
try 
    cd(datdir);
    clear listOftifs;
		clear dir;
    listOftifs = dir('*.tif');
    nfilesmax = numel(listOftifs);
    
if nfilesmax>10
    %data = zeros(256,256,nfilesmax);
    cd(curdir);
    clear data;
    for nfiles = 0:nfilesmax-1
         I = imread([datdir listOftifs(nfiles+1).name]);
         data(:,:,nfiles+1) = I;
        fprintf('Loading file %s \n',[datdir listOftifs(nfiles+1).name])
    end
    sname = sprintf('S%04d',scannum);
    
     [p1 p2 p3] = ind2sub(size(data),find(data>=500));
   %condition of max being greater than 150% of n.n.
   %length(p1)
   numpixmod =0;
%    for qq=1:length(p1)
%    %disp('checking for saturated pixels')
%    nn1 = data(p1(qq)+1,p2(qq),p3(qq));
%    nn2 = data(p1(qq)-1,p2(qq),p3(qq));
%    nn3 = data(p1(qq),p2(qq)+1,p3(qq));
%    nn4 = data(p1(qq),p2(qq)-1,p3(qq));
%    mval = data(p1(qq),p2(qq),p3(qq));
%    
%    if mval>100 && mval>=1.5*nn1 && mval>=1.5*nn2 && mval>=1.5*nn3 && mval>=1.5*nn4
%        disp('saturated pixel found with value')
%        mval
%        disp('changing value to')
%        newval = (1/4)*(nn1+nn2+nn3+nn4)
%  %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
%  %      pause
%        data(p1(qq),p2(qq),p3(qq))=newval;
%  %      scannum
%  %      close(gcf)
%  numpixmod = numpixmod+1
%    end
%    end
    
    
    if padme ==1
        cd(curdir)
     data=fft_padjesse(data,[[1 1],1]);   
    end
    
   % if centdat ==1
   % data=center_array(data);
   % end

 
   
   

  
   
  
    
  %  data = double(sqrt(data));
 	%save(sname,'data');
    
 
% else ind= find(scannum==scannum2);
%                         disp('nothing for')
%                         scannum
% 			scannum2(ind)=0;
% 
%             
% end
end
		cd(curdir);

catch me
    me
        
end
end

params.nn = size(data);
params.binning = [1,1];
