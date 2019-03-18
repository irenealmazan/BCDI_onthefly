%correlate different sets of tif stacks and see how similar they are
clear all; clc;

params.version='Matlab phasing version 1.1 - Nov 2013';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining the current directory
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Data preperation
params.data_dir=dir;            %assumes data is in the same dir as this file
                                %otherwise specify the data directory.
                                % at the moment it will create a save dir with a name that is 
                                % generated (see bottom of page) according to phasing params.  it will be
                                % saved in the directory that this phasing file is in.  this can be over-
                                % ridden by simply specifiying another save_dir.  
params.save_dir=dir;

                                % the data files and background files to load.  if no bg files are needed
                                % put back={}.  Accepts .spe files and also matlab save files .mat.  if
                                % matlab is used it is assumed the array is 3d.  all other options work the
                                % same.

params.files={'35_a.tif'};%
params.back={};%if no bg required, leave empty {}.  don't use {''}.

params.itnum = 1;               %change this for each recon of same data set
 
params.seq_let='';             %sequence letter, used in save ouput. can be any string
params.comments='';      %write comments here and this will be saved in the params file

params.binning=[1,1];           %binning for x and y

params.skipping=0 ;             %will bin (=0), will skip points (=1)

params.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=0;            %min data threshold for NOV/JUNE.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.

%params.min_data=200;            %mindata for Feb recons                             
                                
params.schot_th=0;           %secondary threshold applied AFTER binning
params.subtract_dc=0;           %leave this as 0 (unless you know what you are doing)
params.no_center=0;             %no centering of data -> no_center=1
params.no_fft_pad=1;            %no padding for fft -> no_fft_pad=1
params.pad_ptych=0;             %leave this as 0 (unless you know what you are doing)
params.no_hist=1;               %plot the histograms of the data (=0)
params.bg_mult=1;               %multiplication to apply to the bg file (if exp time is different)
params.save_data='YES';


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
params.do_2D=0;                 %will take the central slice and do 2D phase retrieval if =1                                
                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No inputs required here.  Doing the data preperation.                     
[iterates(:,:,:,1) params] = load_MP_data(params);   %load the data


params.files={'42.tif'};%

[iterates(:,:,:,2) params] = load_MP_data(params);

params.files={'43_a.tif'};%

[iterates(:,:,:,3) params] = load_MP_data(params);

params.files={'45_a.tif'};%

[iterates(:,:,:,4) params] = load_MP_data(params);

params.files={'49_a.tif'};%

[iterates(:,:,:,5) params] = load_MP_data(params);

params.files={'51_a.tif'};%

[iterates(:,:,:,6) params] = load_MP_data(params);
%% now check for saturated pixels

for z=1:3
	data = iterates(:,:,:,z);

[p1 p2 p3] = ind2sub(size(data),find(data>=1000));
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
 %      pause
       data(p1(qq),p2(qq),p3(qq))=newval;
 %      scannum
 %      close(gcf)
   end
	 end
	 
iterates(:,:,:,z) = data;	 
	 
end

%% max intensities are centered so I can just compare the central slice and it should be lined up

for j=1:3

	ref_four=iterates(:,:,:,j);
	
	figure; imagesc(log(ref_four(:,:,20)));
	
% [p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));
% 
% p1s(j)=p1;
% p2s(j)=p2;
% p3s(j)=p3;

%I think this is probably the best way to do it
corrcoef(iterates(:,:,:,1),iterates(:,:,:,j))

end

%% what does the autocorrelation size look like?

for j=1:3

	ref_four=iterates(:,:,:,j);
	
	figure; imagesc(log(xcorr2(ref_four(:,:,20))));
	
	slice = mean(xcorr2(ref_four(:,:,20)));
	
% [p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));
% 
% p1s(j)=p1;
% p2s(j)=p2;
% p3s(j)=p3;

%I think this is probably the best way to do it
%corrcoef(iterates(:,:,:,1),iterates(:,:,:,j))
x=1:511;
f = fit(x.',slice.','gauss1');
sigs(j) = f.c1;


end
