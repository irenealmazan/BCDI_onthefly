% this scripts test whether FT_Irene and InvFT_Irene are equivalent to fftn
% and ifftn respectively

%onthefly_GA_v2; % creates a rhobest

global X Y Z

data = zeros(size(X));

thscantoplot = zeros(1,size(X,3));

thscantoplot_nosifht = zeros(1,size(X,3));

rock_test = zeros(1,size(X,3));

rock_true = zeros(1,size(X,3));

 dqz_shift_toplot = zeros(1,size(X,3));

%rhobest = NW;
%%{
 for ii = 1:numel(data_exp)
     
    data(:,:,ii) = data_exp(ii).I;

    %inibeam(ii).A = probe;
    %inibeam(ii).A = circshift(probe, round([shift_guess(ii,2)  shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

    %%% THETA POSITIONS
    
    %uncorrect theta positions
    %{
    dth_initials(ii) = data_exp(ii).dth; % it is already distorted
    data_exp(ii).dth_new = data_exp(ii).dth;% + data_exp(ii).dth_delta; 
    data_exp(ii).dqshift = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz]; 
    %}
    
    %correct theta positions
    %{
     data_exp(ii).dqshift =  data_exp(ii).dqshift_delta; 
    %}
    
    thscantoplot(ii) = data_exp(ii).dth + data_exp(ii).dth_delta;
    
    thscantoplot_nosifht(ii) = data_exp(ii).dth;% - data_exp(ii).dth_delta;
    
    dqz_shift_toplot(ii) = data_exp(ii).dqshift(3);
    end
%}

% pixel size on the z direction of the reciprocal space:
qz_pixel_size = abs((dqz_shift_toplot(end) - dqz_shift_toplot(1))/numel(dqz_shift_toplot));

% pixel size on the delta direction
%th_pixel_size = 2*pi/abs((dqz_shift_toplot(end) - dqz_shift_toplot(1)));%(thscantoplot(end)-thscantoplot(1))/numel(thscantoplot);  

 % set up experimental details
%thBragg2 = 34; % in deg, taken from labbook page 19
[Npix,Npiy,depth] = size(NW);
%d2_bragg = camdist * params.lam/(Npix*params.det_px);



% [X,Y,Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
%                     [-Npiy/2:Npiy/2-1]*d2_bragg,...
%                     [-depth/2:depth/2-1].*th_pixel_size);

 % fft space:
[X_recip,Y_recip,Z_recip] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
    [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg),...
    [-depth/2:depth/2-1].*2*qz_pixel_size);
                
% FT_Irene space:

X_recip_Irene = X_recip;

for jj = 1:numel(data_exp)
    X_recip_Irene(:,:,jj) = X_recip(:,:,jj) - data_exp(jj).dqshift(1);%repmat(data_exp(jj).dqshift(1),size(X_recip,1),size(X_recip,2)); 
end

Y_recip_Irene = Y_recip;
Z_recip_Irene = Z_recip;



% conventional fftn/ifftn

%%{
%Psi2_true = fftshift(fftn(rhobest));

Psi2_true = fftshift(fftn(NW));

%Psi2_test = data;%fftshift(fftn(rhobest));

%Pmrho_test = ifftn(fftshift(Psi2_test));%fftshift(ifftn(fftshift(Psi2_test)));%


%  test that we obtain the same object with the ifftn(fftn()) function
%%{
%figure; 

for jj=1:size(Psi2_true,3)
    
    %thscantoplot(jj) = data_exp(jj).dth;
    
    %Psi2_test_mod(:,:,jj) = sqrt(Psi2_test(:,:,jj).*conj(Psi2_test(:,:,jj)));
    
    %Psi2_test_mod(:,:,jj) = sqrt(Psi2(:,:,jj).*conj(Psi2(:,:,jj)));
    
    %rock_test(jj) = sum(sum(Psi2_test_mod(:,:,jj),1));
    
    Psi2_true_mod(:,:,jj) = sqrt(Psi2_true(:,:,jj).*conj(Psi2_true(:,:,jj)));
    
    rock_true(jj) = sum(sum(Psi2_true_mod(:,:,jj),1));
    
    %{
    subplot(221);
    imagesc(abs(rhobest(:,:,jj)));
    colorbar;
    axis image;
    title(['rhobest ' num2str(jj)]);
    
    subplot(222);imagesc(angle(rhobest(:,:,jj)));
    colorbar;
    axis image;
    title(['rhobest phase val ' num2str(jj)]);

    subplot(223);
    imagesc(abs(Pmrho_test(:,:,jj)));
    colorbar;
    axis image;
    title(['conventional IFT( FT (rhobest) ) absolute val ' num2str(jj)]);
    
    subplot(224);
    imagesc(angle(Pmrho_test(:,:,jj)));
    colorbar;
    axis image;
    title(['conventional IFT( FT (rhobest) ) phase val ' num2str(jj)]); 
    
    pause(.5); 
    %}

end

%}




% Irene's FT and IFT

%%{
[Psi2,Psi2_mod,rock_Psi2] = FT_Irene(data,data_exp,NW);


%[Pmrho] = InvFT_Irene(data,data_exp,Psi2_mod); % only a test to see what happens with a purely real bragg peak
[Pmrho] = InvFT_Irene(data,data_exp,Psi2);
%}
% comparison of the FT
 
 %{
 figure;
 
 for jj=28%1:size(Psi2,3) 
    
    
    subplot(221);
    imagesc(abs(Psi2_test(:,:,jj)));
    %imagesc(X_recip(1,:,1),Y_recip(:,1,1),abs(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 fftn(rhobest) ' num2str(jj)]);
    
    subplot(222);imagesc(angle(Psi2_test(:,:,jj)));
    %subplot(222);imagesc(X_recip(1,:,1),Y_recip(:,1,1),angle(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 fftn(rhobest) phase val ' num2str(jj)]);

    subplot(223);
    imagesc(abs(Psi2(:,:,jj)));
    %imagesc(X_recip(1,:,1)+data_exp(jj).dqshift(1),Y_recip(:,1,1),abs(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) ) absolute val ' num2str(jj)]);
    
    subplot(224);
    imagesc(angle(Psi2(:,:,jj)));
    %imagesc(X_recip(1,:,1)+data_exp(jj).dqshift(1),Y_recip(:,1,1),angle(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    pause(.5); 

 end

%}


%%{
figure;
 
 for jj=125%1:size(Psi2,1)%1:size(Psi2,3) 
    
    subplot(221);
    imagesc((abs(squeeze(Psi2_true(jj,:,:)))));
    %imagesc(X_recip(1,:,1),Y_recip(:,1,1),abs(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 true ' num2str(jj)]);
    
%     subplot(232);
%     imagesc((abs(squeeze(Psi2_test(jj,:,:)))));
%     %imagesc(X_recip(1,:,1),Y_recip(:,1,1),abs(Psi2_test(:,:,jj)));
%     colorbar;
%     axis image;
%     title(['Psi2 fftn(rhobest) ' num2str(jj)]);
    
    subplot(222);
    imagesc((abs(squeeze(Psi2(jj,:,:)))));
    %imagesc(X_recip(1,:,1)+data_exp(jj).dqshift(1),Y_recip(:,1,1),angle(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) )  ' num2str(jj)]); 
    
    subplot(223);imagesc(angle(squeeze(Psi2_true(jj,:,:))));
    %subplot(222);imagesc(X_recip(1,:,1),Y_recip(:,1,1),angle(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 true phase val ' num2str(jj)]);

    
    
    
%     subplot(235);imagesc(angle(squeeze(Psi2_test(jj,:,:))));
%     %subplot(222);imagesc(X_recip(1,:,1),Y_recip(:,1,1),angle(Psi2_test(:,:,jj)));
%     colorbar;
%     axis image;
%     title(['Psi2 fftn(rhobest) phase val ' num2str(jj)]);

  
    
    subplot(224);
    imagesc(angle(squeeze(Psi2(jj,:,:))));
    %imagesc(X_recip(1,:,1)+data_exp(jj).dqshift(1),Y_recip(:,1,1),angle(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    pause(.5); 

 end
%}

% comparison between IFT_Irene(FT_Irene()) and ifftn(fftn())
%{
figure; 

for jj=100%80:120 
    
    
    
    subplot(231);
    imagesc(abs(squeeze(rhobest(jj,:,:))));
    colorbar;
    axis image;
    title(['rhobest ' num2str(jj)]);
   
    
    subplot(232);
    imagesc(abs(squeeze(Pmrho_test(jj,:,:))));
    colorbar;
    axis image;
    title(['Conventional IFT(data ) absolute val ' num2str(jj)]);
    
    
    subplot(233);
    imagesc(abs(squeeze(Pmrho(jj,:,:))));
    colorbar;
    axis image;
    title(['Irene IFT( FT (rhobest) ) absolute val ' num2str(jj)]);
    
    subplot(234);imagesc(angle(squeeze(rhobest(jj,:,:))));
    colorbar;
    axis image;
    title(['rhobest phase val ' num2str(jj)]);
   
     subplot(235);
    imagesc(angle(squeeze(Pmrho_test(jj,:,:))));
    colorbar;
    axis image;
    title(['Irene IFT( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    subplot(236);
    imagesc(angle(squeeze(Pmrho(jj,:,:))));
    colorbar;
    axis image;
    title(['Irene IFT( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    pause(.5); 

end
%}
figure;

jj= 25;%round(size(rhobest,3)/2);

 subplot(221);
 imagesc(abs(rhobest(:,:,jj)));
 colorbar;
 axis image;
 title(['rhobest ' num2str(jj)]);
 
 subplot(222);imagesc(angle(rhobest(:,:,jj)));
 colorbar;
 axis image;
 title(['rhobest phase val ' num2str(jj)]);
 
 subplot(223);
 imagesc(abs(Pmrho(:,:,jj)));
 colorbar;
 axis image;
 title(['Irene IFT( FT (rhobest) ) absolute val ' num2str(jj)]);
 
 subplot(224);
 imagesc(angle(Pmrho(:,:,jj)));
 colorbar;
 axis image;
 title(['Irene IFT( FT (rhobest) ) rhobest phase val ' num2str(jj)]);
 
 figure;
 plot(thscantoplot_nosifht,rock_true,'or');
 hold on;
 plot(thscantoplot,rock_Psi2,'*k'); 
 %plot(thscantoplot,(rock_test),'ob');
 legend('true rock curve','Irene rock curve');
 
%  figure;
%  plot(dqz_shift_toplot,rock_Psi2,'*r'); 
%  hold on; 
%  plot(dqz_shift_toplot,(rock_test),'ob');
%  
 
 
 figure;
 subplot(121);
 di(Psi2_true,-0.01,'g',X_recip,Y_recip,Z_recip);
 axis image
 title('FT true FT(NW)');
 
%  subplot(132);
%  di(Psi2_test,-0.01,'b',X_recip,Y_recip,Z_recip);
%  axis image
%  title('FT conventional');
 
 subplot(122);
 di(Psi2,-0.01,'r',X_recip_Irene,Y_recip_Irene,Z_recip_Irene);
 axis image;
 title('FT Irene');
 
 
 
 figure;
 subplot(121);
 di(NW,-0.01,'g',X,Y,Z);
 axis image;
 title('True object');
 
 subplot(122);
 di(rhobest,-0.09,'r',X,Y,Z);
 axis image;
 title('IFT(FT) Irene');
 
%  subplot(122);
%  di(Pmrho,-0.09,'r',X,Y,Z);
%  axis image;
%  title('IFT(FT) Irene');
 
%  subplot(133);
%  di(Pmrho_test,-0.01,'r',X,Y,Z);
%  axis image;
%  title('Conventional IFT(FT)');
%  
%  
 
return;

%%% visual representation
 %{
 figure;
 
 jj = 1;  
 
Rvect = [0 0 1];
Rcenter = [0 0 0];
 
% plane 1 of the beam
[ grid(jj).yy1,  grid(jj).xx1 ,  grid(jj).zz1] = meshgrid([-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                    [-Npix/2:Npix/2-1].*2*pi/(Npix*d2_bragg), [0]);

% sets the cartesian coordinates 
grid(jj).points2d_1 = [  grid(jj).xx1(:)  grid(jj).yy1(:)];


% create a 3rd dimension and convert the pixel number in inverse angstroms:
grid(jj).points3d_1 = [  grid(jj).points2d_1 grid(jj).zz1(:)]; 


% reshape the mesh grid which is now centered around the beam maximum
% intensity and in microns units:
grid(jj).xx1_microns = reshape(grid(jj).points3d_1(:,1),size(grid(jj).xx1));
grid(jj).yy1_microns = reshape(grid(jj).points3d_1(:,2),size(grid(jj).xx1));
grid(jj).zz1_microns = reshape(grid(jj).points3d_1(:,3),size(grid(jj).xx1));

% plot
hold on;
h(jj) = surf(grid(jj).xx1_microns, grid(jj).yy1_microns, grid(jj).zz1_microns);
h(jj).EdgeColor = 'none'; % set properties of the plot
%h(jj).CData = abs(beam); % use the intensity of the beam to color the surface
axis image

% rotation
rotate(h(jj), Rvect, thBragg, Rcenter);

hold off;

% extract the rotated coordinates in the lab frame:

rotpos.x = h(jj).XData;
rotpos.y = h(jj).YData;
rotpos.z = h(jj).ZData;


[X_recip_toplot Y_recip_toplot Z_recip_toplot] = meshgrid([-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                    [-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg),...
                    [-depth/2 depth/2-1].*2*qz_pixel_size);
                
  
               subplot(131)
hold on;
h_Psi2 =  di(Psi2_test,-0.01,'g',X_recip,Y_recip,Z_recip);

[X_recip_toplot Y_recip_toplot Z_recip_toplot] = meshgrid([-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg), ...
                    [-Npix/2 Npix/2-1].*2*pi/(Npix*d2_bragg),...
                    [-depth/2 depth/2-1].*2*pi/(depth*d2_bragg));

scatter3(X_recip_toplot(:),Y_recip_toplot(:),Z_recip_toplot(:));
 
                
for jj=1:numel(data_exp)
    
  
    
    Psi2_toplot(:,:,jj) = temp;
   %translate the detector plane from dq_shift:
   trrotpos(jj).x = rotpos.x + data_exp(jj).dqshift(1);
   trrotpos(jj).y = rotpos.y + data_exp(jj).dqshift(2);
   trrotpos(jj).z = rotpos.z + data_exp(jj).dqshift(3);    
   
   
   % plot
  subplot(131);
  hold on;
   h(jj) = surf(trrotpos(jj).x, trrotpos(jj).y, trrotpos(jj).z);
   h(jj).EdgeColor = 'none'; % set properties of the plot
   h(jj).CData = abs(Psi2_mod(:,:,jj)); % use the intensity of the beam to color the surface
   alpha(h(jj),0.3);
   
   drawnow;

   axis image;
   %hold off;
   
   subplot(132);
   imagesc(Psi2_test_mod(:,:,jj));
   axis image;
   title('Conventional FT');
   colorbar;
   
   subplot(133);
   imagesc(Psi2_mod(:,:,jj));
   axis image;
   title('Irene FT');
   colorbar;
   
end
%}

 