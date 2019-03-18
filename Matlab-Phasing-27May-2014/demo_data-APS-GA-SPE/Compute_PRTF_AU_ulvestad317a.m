%load up the data
clear all; close all; clc;
load('/Users/b228289/Box Sync/PdH grains/Ulvestad317adsets/S0204.mat')

%%
%play_data(data)
data = sqrt(data);
%%
%now with data, can do prtf
amps = amps/max(amps(:));
amps(amps<0.2)=0;
ref = amps.*exp(1i.*phs);

%no pcdi
ref_four = fftshift(fftn(ref));


[p1 p2 p3] = ind2sub(size(ref_four),find(abs(ref_four)==max(max(max(abs(ref_four))))));

[p11 p21 p31] = ind2sub(size(data),find(data==max(max(max(data)))));

%%
for zz=p3-2:1:p3+2
figure;
zoom(2);
imagesc(log(abs(ref_four(:,:,zz))));
end
%%
for zz=p31-2:1:p31+2
figure;
zoom(2);
imagesc(log(flipud(data(:,:,zz))));
end
%%
figure; imagesc(log(abs(ref_four(:,:,p3)))); zoom(3)
figure; imagesc(log(data(:,:,p31))); zoom(3)
figure; imagesc(log(flipud(data(:,:,p31)))); zoom(3);
%figure; imagesc(log(fliplr(data(:,:,p31)))); zoom(3);
%%
%now choose the slices, check for being flip up down
rcnst_abs=(abs(ref_four(:,:,p3)));
meas_abs=flipud(data(:,:,p31));
%meas_abs=flipud(datan(:,:,24));

%rescale
rcnst_abs = rcnst_abs/max(rcnst_abs(:));
meas_abs = meas_abs/max(meas_abs(:));

%check to see that they look similar
figure; imagesc(log10(rcnst_abs));
figure; imagesc(log10(meas_abs));

%% make sure both are centered
rcnst_abs(end/2,end/2)
meas_abs(end/2,end/2)

rcnst_abs=center_array_data(rcnst_abs);
rcnst_abs(end/2,end/2)
%% may need to align as best you can
% [h k l] = register_3d_reconstruction(meas_abs,rcnst_abs);
% next_c = sub_pixel_shift(rcnst_abs,h,k,l);
% 
% figure; imagesc(log10(rcnst_abs));
% figure; imagesc(log10(abs(next_c)));
% 
% rcnst_abs = abs(next_c);
%% Do prtf for 5 slices surrounding the max in the reconstructed set

%rcnst_abs=abs(ref_four(:,:,p3));

for hr=0;%-2:1:2

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
  
 % imagesc(temp22);
 % pause(0.2)

  %figure(858); imagesc(log10(1+abs( temp22 .* rcnst_abs )))
  
 % prtf(kk) = sum(sum(temp22 .* rcnst_abs)) / (1E-9 + sum(sum(temp22 .* meas_abs)));
  
 prtf(kk) = sum(sum(temp22 .* meas_abs)) / (1E-9 + sum(sum(temp22 .* rcnst_abs)));
  
  kk = kk+1;

end
toc

prtf2(:,hr+3)=prtf';

end
%%
figure; plot(prtf)

%%
figure; plot(prtf2(:,3))
%%
cd ~/Documents/Ulvestad317/;
[specscan, errors] = openspec('Ulvestad317a.spec', 204);
camdist = specscan.motor_positions(27)/1000; % in meters
%convert it to Q inverse angstrom
maxpix = 128;
max2th = atand((1:maxpix)*(55e-6)/camdist); %pixel size/cam dist
angres = (2*pi)./((4*pi)./(1.377)*sind(max2th)); %in angstrom
maxQ = ((4*pi)./(1.377)*sind(max2th));

cd ~/Desktop/;
%smooth the prtf before plotting
figure; plot(maxQ(1:64),smooth(prtf(1:64),5));
xlabel('Max Q (A^-1)')
ylabel('PRTF')
%% keep it in pixels then use angres 
figure; plot((1:64),smooth(prtf(1:64),5));


%% in angstrom resolution
figure; plot(angres(1:64),smooth(prtf(1:64),5));
xlabel('Angstrom resolution')
ylabel('PRTF')
set(gcf, 'xdir','reverse')
%saveas(1,'prtf.eps','psc2')
%%
figure; plot([1:length(prtf(1:100))]*(1.3047e-4),smooth(prtf(1:100),1));
%--------------------------------------------------------------------------
figure; semilogx(prtf,'-or','LineWidth',2,'MarkerSize',4); 
xlim([1 length(prtf)])
%ylim([1E-7 1E5])
%--------------------------------------------------------------------------
clear('xx','yy','annuls_smplng','jj_max','loop_range','jj','kk')
clear('temp*','o_diam_circ','i_diam_circ')
%--------------------------------------------------------------------------

