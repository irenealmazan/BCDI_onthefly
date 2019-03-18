%load up the data

params2.version='Matlab phasing version 1.1 - Nov 2013';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determining the current directory
name_of_this_file='Matlabphasing_ver1_1';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
params2.this_dir=dir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Data preperation
params2.data_dir=dir;            %assumes data is in the same dir as this file
                                %otherwise specify the data directory.
                                % at the moment it will create a save dir with a name that is 
                                % generated (see bottom of page) according to phasing params.  it will be
                                % saved in the directory that this phasing file is in.  this can be over-
                                % ridden by simply specifiying another save_dir.  
params2.save_dir=dir;

                                % the data files and background files to load.  if no bg files are needed
                                % put back={}.  Accepts .spe files and also matlab save files .mat.  if
                                % matlab is used it is assumed the array is 3d.  all other options work the
                                % same.

params2.files={'1901.tif'};%
params2.back={};%if no bg required, leave empty {}.  don't use {''}.

params2.itnum = 1;               %change this for each recon of same data set
 
params2.seq_let='';             %sequence letter, used in save ouput. can be any string
params2.comments='';      %write comments here and this will be saved in the params file

params2.binning=[1,1];           %binning for x and y

params2.skipping=0 ;             %will bin (=0), will skip points (=1)

params2.aliens=[0];           
                                %set aliens=0 for no removal, otherwise input them as    
                                %aliens=[[x0,y0,z0,x1,y1,z1],[x2,y2,z2,x3,y3,z3]]
                                %will remove two instances of aliens given by the pairs
                                % #0 and #1 and another given by #2,#3. accepts as
                                % many as you like.  Points are the same as winview
                         
params2.min_data=3;            %min data threshold for NOV/JUNE.  the 
                                %threshold is applied to each data set BEFORE
                                %addition and binning.
                      
                                
params2.schot_th=0;           %secondary threshold applied AFTER binning
params2.subtract_dc=0;           %leave this as 0 (unless you know what you are doing)
params2.no_center=0;             %no centering of data -> no_center=1
params2.no_fft_pad=0;            %no padding for fft -> no_fft_pad=1
params2.pad_ptych=0;             %leave this as 0 (unless you know what you are doing)
params2.no_hist=1;               %plot the histograms of the data (=0)
params2.bg_mult=1;               %multiplication to apply to the bg file (if exp time is different)
params2.save_data='YES';


params2.data_only='NO';        
                                %will save or exit after data prep. ='NO' no save and no exit,         
                                %'YES' will save and continue, 'YES-EXIT'
                                %will save and then stop the script

params2.nnc=[0,0,0,0,15,15];
                                % initial cropping of data before binning.
                                %eg. nnc=[0,0,-10,-10,5,5] will do nothing to x,
                                %will crop 10 pixels off each end in y and will pad
                                %5 pixels to each end in z.  set nnc=0 to do
                                %nothing.
params2.do_2D=0;                 %will take the central slice and do 2D phase retrieval if =1                                
                                
                   
[data params2] = load_MP_data(params2);   %load the data


%%
%now with data, can do prtf


%make sure to form rcnst_abs by taking FT(amp>0.2 * exp[i phase*amp>0.2))
%load in the data and the params that has the reconstruction
iterates = params.pnm_avg;

ref = iterates(:,:,:,9);
datan = data;

%ref_amp=abs(ref).*(abs(ref)>0.001);
%ref_ph = angle(ref).*(abs(ref)>0.001);

%ref_four = fftshift(fftn(ref_amp.*exp(1i.*ref_ph)));
ref_four = fftshift(fftn(ref));

%find where the maxes are

[p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));

[p11 p21 p31] = ind2sub(size(data),find(data==max(max(max(data)))));


for zz=p3-2:1:p3+2
figure;
imagesc(abs(ref_four(:,:,zz)));
end

for zz=p31-2:1:p31+2
figure;
imagesc(flipud(data(:,:,zz)));
end

figure; imagesc(abs(ref_four(:,:,35))); zoom(3)
figure; imagesc(flipud(data(:,:,34))); zoom(3)


%now choose the slices
rcnst_abs=abs(ref_four(:,:,50));
meas_abs=flipud(data(:,:,44));
%meas_abs=flipud(datan(:,:,24));

%rescale
rcnst_abs = rcnst_abs./max(max(rcnst_abs));
meas_abs = meas_abs./max(max(meas_abs));

%check to see that they look similar
figure; imagesc(rcnst_abs);
figure; imagesc(meas_abs);

%% Do prtf for 5 slices surrounding the max in the reconstructed set

rcnst_abs=abs(ref_four(:,:,p3));

for hr=-2:1:2

meas_abs=flipud(data(:,:,p31+hr));

rcnst_abs = rcnst_abs./max(max(rcnst_abs));
meas_abs = meas_abs./max(max(meas_abs));

[Nx Ny] = size(rcnst_abs); %should only be 2D

%#ok<*BDSCI>
%--------------------------------------------------------------------------
%keep this for the time being, compare speeds for defining it here vs passing it
[xx,yy] = meshgrid(1:Nx,1:Ny);
%--------------------------------------------------------------------------
%the q = 0 value:
jj = 1; kk = 1;

o_diam_circ = (((Nx/2+1-xx)./jj).^2 + ((Ny/2+1-yy)./jj).^2 < 1);

prtf(kk) = sum(sum(o_diam_circ .* rcnst_abs)) / ...
                      (1E-9 + sum(sum(o_diam_circ .* meas_abs)));
%--------------------------------------------------------------------------
jj_max = round(Nx/2); 
annuls_smplng = 1;

%after q = 0, which corresponds to jj = 1, we now need to loop from:
loop_range = ((1 + annuls_smplng) : annuls_smplng : jj_max);
%--------------------------------------------------------------------------
prtf(2 : (1 + length(loop_range))) = 0; 
%--------------------------------------------------------------------------
tic
kk = 2;
for jj = loop_range
    
  o_diam_circ = (((Nx/2+1-xx)./jj).^2 + ((Ny/2+1-yy)./jj).^2 < 1);
  i_diam_circ = (((Nx/2+1-xx)./(jj-1)).^2 + ((Ny/2+1-yy)./(jj-1)).^2 < 1);

%i_diam_circ = (((Nx/2+1-xx)./(jj-5)).^2 + ((Ny/2+1-yy)./(jj-5)).^2 < 1);

  temp22 = o_diam_circ - i_diam_circ;
  
%   imagesc(temp22);
%   pause

  %figure(858); imagesc(log10(1+abs( temp22 .* rcnst_abs )))
  
%   prtf(kk) = sum(sum(temp22 .* rcnst_abs)) / (1E-9 + sum(sum(temp22 .* meas_abs)));
%   
 prtf(kk) = sum(sum(temp22 .* meas_abs)) / (1E-9 + sum(sum(temp22 .* rcnst_abs)));
  
  kk = kk+1;

end
toc

prtf2(:,hr+3)=prtf';

end

hold on
for g=1:5
	plot(prtf2(1:35,g))
	pause
end

%find the one that looks the best then check to make sure it actually looks
%similar
figure; imagesc(rcnst_abs);
meas_abs=flipud(data(:,:,p31));
meas_abs = meas_abs./max(max(meas_abs));
figure; imagesc(meas_abs);

%smooth the prtf before plotting
figure; plot([1:90]*(1.3047e-4),smooth(prtf(1:90),3));
xlabel('Q (Ansgtrom^{-1})')
ylabel('PRTF')

saveas(1,'prtf.eps','psc2')

figure; plot([1:length(prtf(1:30))]*(1.3047e-4),smooth(prtf(1:30),1));
%--------------------------------------------------------------------------
figure; semilogx(prtf,'-or','LineWidth',2,'MarkerSize',4); 
xlim([1 length(prtf)])
%ylim([1E-7 1E5])
%--------------------------------------------------------------------------
clear('xx','yy','annuls_smplng','jj_max','loop_range','jj','kk')
clear('temp*','o_diam_circ','i_diam_circ')
%--------------------------------------------------------------------------

