%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlabphasing_ver1.1.m 
%2D/3D phase retrieval - Jesse Clark, LCN, UCL October-November 2010
%                       jesse.clark@ucl.ac.uk, jessenclark@gmail.com
clear all; clc;
params.version='Matlab phasing version 1.1 - Nov 2010';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining the current directory
name_of_this_file='Data_prep_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Data preperation
%params.data_dir='/Users/a4894z/Dropbox/Matlab_phasing_AUG13/GA/';
params.data_dir='/Users/a4894z/Desktop/';%assumes data is in the same dir as this file
                                %otherwise specify the data directory.
                                % at the moment it will create a save dir with a name that is 
                                % generated (see bottom of page) according to phasing params.  it will be
                                % saved in the directory that this phasing file is in.  this can be over-
                                % ridden by simply specifiying another save_dir.  
params.save_dir='/Users/a4894z/Desktop/';

                                % the data files and background files to load.  if no bg files are needed
                                % put back={}.  Accepts .spe files and also matlab save files .mat.  if
                                % matlab is used it is assumed the array is 3d.  all other options work the
                                % same.

params.files={'69.tif','70.tif'};%,'42.tif','43_a.tif','45_b.tif','49.tif','51_a.tif','53_a.tif','55_a.tif','57.tif','58_a.tif'};
                                %specify the files.  they should be one for
                                %each data set.  to combine multiple for a
                                %single data set, use Data_prep_ver1_1.m or
                                %something similiar.
                                
params.back={};                 %one bg file for each data set.  leave blank if already done


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
%% build in data correction for saturated pixels

for ze = 1:params.ntime
	
	data = params.data(:,:,:,ze);

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
	 
	 params.data(:,:,:,ze) = data;
end

%% max intensities are centered so I can just compare the central slice and it should be lined up

initialize_plot_PRL_old;

iterates = params.data;

[p1 p2 p3] = ind2sub(size(iterates(:,:,:,1)),find(abs(iterates(:,:,:,1))==max(max(max(abs(iterates(:,:,:,1)))))));


for j=1:params.ntime

	ref_four=iterates(:,:,:,j);
	
	figure; imagesc(log(ref_four(:,:,p3))); title(params.files(j));
	
% [p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));
% 
% p1s(j)=p1;
% p2s(j)=p2;
% p3s(j)=p3;

%I think this is probably the best way to do it
if j<params.ntime
[r,p]=corrcoef(iterates(:,:,:,1),iterates(:,:,:,j+1));
coefs(j) = r(1,2);
end



end

figure; plot(coefs,'.-');
%figure; plot(latc,'.-');

scannum_new = scannum(coefs>0.7);


%%

data = 0.5*(params.data(:,:,:,1)+params.data(:,:,:,2));

mat2tif(data,[num2str(69) '_a.tif'])

%%

mat2tif(data(:,:,1:81),[num2str(1321) '.tif'])

mat2tif(data(:,:,82:162),[num2str(1322) '.tif'])

mat2tif(data(:,:,163:end),[num2str(1323) '.tif'])

%mat2tif(data,[num2str(45) '_b' '.tif'])


%% save the data
%get the number for the file name
nfiles=numel(params.files);
sname=[];
for qq=1:nfiles

    a=cell2mat(params.files(qq));
    numb=num2str(sscanf(a(strfind(a,'-')+1:numel(a)),'%i'));

    if numel(numb) == 0,numb=num2str(sscanf(a(strfind(a,'_')+1:numel(a)),'%i'));end

    numbs(qq)=str2num(numb);
    
    szn=size(numb);
    if szn(1) > 1,numb=numb(end,:);end
    
    if qq ~= nfiles, sname=[sname,numb,'-'];else sname=[sname,numb];end

    
end
    
mat2tif(data,[num2str(190) '.tif'])


savevtk2scalar(data,[dir,'/',params.seq_data,sname,'.vtk'])