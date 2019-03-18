%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matlabphasing_ver1.1.m 
%2D/3D phase retrieval - Jesse Clark, LCN, UCL October-November 2010
%                       jesse.clark@ucl.ac.uk, jessenclark@gmail.com
clc
clear all
params.version='Matlab phasing version 1.1 - Nov 2013';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining the current directory
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data preperation
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

params.files={'dsring.mat'};%
params.back={'dsbkg.mat'};%if no bg required, leave empty {}.  don't use {''}.

params.itnum = 1;               %change this for each recon of same data set
 
params.seq_let='DSring';             %sequence letter, used in save ouput. can be any string
params.comments='';      %write comments here and this will be saved in the params file

params.binning=[1,1];           %binning for x and y

params.skipping=0 ;             %will bin (=0), will skip points (=1)

params.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params.min_data=200;            %min data threshold.  below this is removed.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.
                                
params.schot_th=0;           %secondary threshold applied AFTER binning
params.subtract_dc=0;           %leave this as 0 (unless you know what you are doing)
params.no_center=1;             %no centering of data -> no_center=1
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
%% No inputs required here.  Doing the data preperation.                     
[data params] = load_MP_data(params);   %load the data

for i=1:566
    imagesc(data(:,:,i))
    pause
    i
end

%write this into a movie

% Prepare the new file.
    vidObj = VideoWriter('dsring.avi');
    open(vidObj);
 
    % Create an animation.

 
 for i=1:566
        imagesc(data(:,:,i))
       % Write each frame to the file.
       currFrame = getframe;
       writeVideo(vidObj,currFrame);
 end
  
    % Close the file.
    close(vidObj);


%this is some guys code and works well

% http://physics.georgetown.edu/matlab/tutorial.html

 %vidObj = VideoWriter('dsring_filt.avi');
 %open(vidObj);

 
%will need to do up to phase transition, the phase transition, then back
%require different parameters for max jumps, etc, how to identify

%up to phase transition is frame 185
 
clear pos;
rnew=0;
hold on;
for j=1:185 %do up through the phase transition
b = bpass(data(:,:,j),1,8); %filter the image
colormap('gray'),image(b);
pk = pkfnd(b,100,10); %find the peaks, min brightness, then size
%min bright of 40 finds a lot
%min bright of 100 only finds 8
if isempty(pk)==0
cnt = cntrd(b,pk,10); %fit centroid of size 10, last option is to show centroid calc
 
colormap('gray'),image(b);
plot(cnt(:,1),cnt(:,2),'rx','Markersize',18)
%pause
 %currFrame = getframe;
 %      writeVideo(vidObj,currFrame);

%position list should be 
%   (x)      (y)      (t)
%   ;     pos = 3.60000      5.00000      0.00000
%   ;           15.1000      22.6000      0.00000
%   ;           4.10000      5.50000      1.00000 
%   ;           15.9000      20.7000      2.00000
%   ;           6.20000      4.30000      2.00000


%this will work for the first time. Then need to bump down to empty rows to
%start next
pos(1+rnew:length(cnt(:,1))+rnew,1:2) = cnt(:,1:2);

pos(1+rnew:length(cnt(:,1))+rnew,3) = j; 

[rnew cnew] = size(pos);

end


end

%now that I have the position list, use the track function
maxdisp=3;
param.mem=5; %once particle is lost, count it as a new particle after 3 steps
param.good=10; %eliminate with fewer than 3 valid positions
param.dim=2;
param.quiet=0;

res = track(pos,maxdisp,param);

%res is in format x y t id#

%compare particle by particle, find params that give particles
%that don't move in the vertical
for partid = 1:res(end,4)

ind = find(res(:,4) == partid);

%plot the particle
%plot(res(ind,1),res(ind,2),'rx')

%determine its velocity
vel = gradient(res(ind,2),res(ind,3));

vbar(partid) = mean(vel);
vsig(partid) = std(vel);

nlocs(partid) = length(ind);%succesive locations


%pause; clf;

end

%close(vidObj);

figure;
hist(nlocs,100)
xlabel('number of frames tracked')
ylabel('number of particles')

mean(nlocs)
std(nlocs)

mean(vbar)
mean(vsig)