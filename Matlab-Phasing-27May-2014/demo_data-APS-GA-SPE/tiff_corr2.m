%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MultiplePhasing_ver2.0.m 
%2D/3D phase retrieval - Jesse Clark,  October 2013
%                       jesclark@stanford.edu, jessenclark@gmail.com
clear all; clc; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining the current directory
name_of_this_file='MultiplePhasing_ver2_0.m';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Data preperation


%snums=[196 197 199 200 203 205 207 209 210 212];

params.save_dir=dir;

                                % the data files and background files to load.  if no bg files are needed
                                % put back={}.  Accepts .spe files and also matlab save files .mat.  if
                                % matlab is used it is assumed the array is 3d.  all other options work the
                                % same.

params.data_dir={'/Users/andrewulvestad/Dropbox/Matlab_phasing_AUG13/GA/'};

                                %data directories for the data.  need to specify for each data set  

%june post cycling exp 4
%params.files={'196.tif','197.tif','199.tif','200.tif','203.tif','205.tif','207.tif','209.tif','210.tif','212.tif'};

%june exp 1 disconnection
%params.files={'35_a.tif','41.tif','42.tif','43_a.tif','45_a.tif','49_a.tif','51_a.tif','53_a.tif','55_a.tif','57.tif','58_a.tif'};
                                %specify the files.  they should be one for
                                %each data set.  to combine multiple for a
                                %single data set, use Data_prep_ver1_1.m or
                                %something similiar.
                    
params.files={'35_a.tif','41.tif','42.tif','43.tif','45.tif','49.tif','51.tif','53.tif','55.tif','57.tif','58.tif'}; %for mult phasing schematic
                                
params.back={};                 %one bg file for each data set.  leave blank if already done


params.seq_let='J1';             %sequence letter, used in save ouput. can be any string
params.comments='june experiment 1 disconnection';      %write comments here and this will be saved in the params file

params.binning=[1,1];           %binning for x and y

params.skipping=0 ;             %will bin (=0), will skip points (=1)

params.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=1;            %min data threshold.  below this is removed.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.

																%always center june data                                
params.no_center=0;             %no centering of data -> no_center=1
params.no_fft_pad=1;            %no padding for fft -> no_fft_pad=1



params.data_only='NO';        
                                %will save or exit after data prep. ='NO' no save and no exit,         
                                %'YES' will save and continue, 'YES-EXIT'
                                %will save and then stop the script

params.nnc=[0,0,0,0,0,0];
                                % initial cropping of data before binning.
                                %eg. nnc=[0,0,-10,-10,5,5] will do nothing to x,
                                %will crop 10 pixels off each end in y and will pad
                                %5 pixels to each end in z.  set nnc=0 to do
                                %nothing.
                               
                
[params] = load_Multi_data(params);   %load the data

%some experimental parameters

%june 1 disconnection
latc = [8.069126991
8.091579851
8.050852155
8.144051458
8.141
8.131157344 %49_a
8.162428088
8.161885343
8.1613
8.163946327
8.175714218];

%june 4 post cycling
latc = [8.173212187
8.179662459
8.175844893
8.169025791
8.114736231
8.098563875
8.143313028
8.174397507
8.193170026
8.200700122];

voltc = [3.57
4.08
4.42
4.68
4.74
4.68
4.6
4.42
3.95
3.66];

%% now check for saturated pixels

for z=1:params.ntime
	data = params.data(:,:,:,z);

[p1 p2 p3] = ind2sub(size(data),find(data>=2000));
   %condition of max being greater than 150% of n.n.
   %length(p1)
   for qq=1:length(p1)
   
   nn1 = data(p1(qq)+1,p2(qq),p3(qq));
   nn2 = data(p1(qq)-1,p2(qq),p3(qq));
   nn3 = data(p1(qq),p2(qq)+1,p3(qq));
   nn4 = data(p1(qq),p2(qq)-1,p3(qq));
   mval = data(p1(qq),p2(qq),p3(qq));
   
   if mval>50 && mval>=1.5*nn1 && mval>=1.5*nn2 && mval>=1.5*nn3 && mval>=1.5*nn4
       disp('saturated pixel found with value')
       mval
       disp('changing value to')
       newval = (1/4)*(nn1+nn2+nn3+nn4)
 %      figure; imagesc(data(:,:,p3(qq))); zoom(3);
 z
       pause
       data(p1(qq),p2(qq),p3(qq))=newval;
 %      scannum
 %      close(gcf)
   end
	 end
	 
params.data(:,:,:,z) = data;	 
	 
end

%z=8;
%data = params.data(:,:,:,z) ;	

%mat2tif(data,[num2str(209) '.tif'])



%% max intensities are centered so I can just compare the central slice and it should be lined up

iterates = params.data;

[p1 p2 p3] = ind2sub(size(iterates(:,:,:,1)),find(abs(iterates(:,:,:,1))==max(max(max(abs(iterates(:,:,:,1)))))));


for j=1:params.ntime

	ref_four=iterates(:,:,:,j);
	
	figure; imagesc(log(ref_four(:,:,p3)));
	
% [p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));
% 
% p1s(j)=p1;
% p2s(j)=p2;
% p3s(j)=p3;

%I think this is probably the best way to do it
if j<params.ntime
[r,p]=corrcoef(iterates(:,:,:,j),iterates(:,:,:,j+1));
coefs(j) = r(1,2);
end



end

%figure; plot(coefs,'.-');
%figure; plot(latc,'.-');

%% what does the autocorrelation size look like?

for j=1:params.ntime

	ref_four=iterates(:,:,:,j);
	
	autco2=abs(fftshift(ifftn(fftshift(ref_four.^2))));
	
	figure; imagesc(autco2(:,:,p3));
	
	slice = mean(autco2(:,:,p3));
	
% [p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));
% 
% p1s(j)=p1;
% p2s(j)=p2;
% p3s(j)=p3;

%I think this is probably the best way to do it
%corrcoef(iterates(:,:,:,1),iterates(:,:,:,j))
x=1:256;
f = fit(x.',slice.','gauss1');
sigs(j) = f.c1;


end

figure; plot(sigs,'.-');
