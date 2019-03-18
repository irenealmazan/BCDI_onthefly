% this scripts test whether FT_Irene and InvFT_Irene are equivalent to fftn
% and ifftn respectively

%onthefly_GA_v2; % creates a rhobest

rhobest = NW;
%%{
 for ii = 1:numel(data_exp)
     
    data(:,:,ii) = data_exp(ii).I;

    %inibeam(ii).A = probe;
    %inibeam(ii).A = circshift(probe, round([shift_guess(ii,2)  shift_guess(ii,1) shift_guess(ii,3)]/d2_bragg));

    %%% THETA POSITIONS
    
    %uncorrect theta positions
    %{
    dth_initials(ii) = data_exp(ii).dth;
    data_exp(ii).dth_new = data_exp(ii).dth; 
    data_exp(ii).dqshift = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz]; 
    %}
    
    %correct theta positions
    %{
     data_exp(ii).dqshift =  data_exp(ii).dqshift_delta; 
    %}
    
    
    end
%}

    
 % set up experimental details
thBragg2 = 34; % in deg, taken from labbook page 19
[Npix,Npiy,depth] = size(rhobest);
d2_bragg = camdist * params.lam/(Npix*params.det_px);



[X,Y,Z] = meshgrid([-Npix/2:Npix/2-1]*d2_bragg, ...
                    [-Npiy/2:Npiy/2-1]*d2_bragg,...
                    [-depth/2:depth/2-1]*d2_bragg);




% conventional fftn/ifftn

%%{
Psi2_test = fftshift(fftn(rhobest));

Pmrho_test = ifftn(fftshift(Psi2_test));

figure; 

for jj=1:size(Pmrho_test,3)
    
    %thscantoplot(jj) = data_exp(jj).dth;
    
    Psi2_test_mod = sqrt(Psi2_test(:,:,jj).*conj(Psi2_test(:,:,jj)));
    
    rock_test(jj) = sum(sum(Psi2_test_mod,1));
    
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

end

%}




% Irene's FT and IFT

[Psi2,Psi2_mod,rock_Psi2] = FT_Irene(data,data_exp,rhobest);

[Pmrho] = InvFT_Irene(data,data_exp,Psi2);

figure; 

for jj=100:size(Pmrho,1) 
    
    %thscantoplot(jj) = data_exp(jj).dth;
    
    subplot(221);
    imagesc(abs(squeeze(rhobest(jj,:,:))));
    colorbar;
    axis image;
    title(['rhobest ' num2str(jj)]);
    
    subplot(222);imagesc(angle(squeeze(rhobest(jj,:,:))));
    colorbar;
    axis image;
    title(['rhobest phase val ' num2str(jj)]);

    subplot(223);
    imagesc(abs(squeeze(Pmrho(jj,:,:))));
    colorbar;
    axis image;
    title(['Irene IFT( FT (rhobest) ) absolute val ' num2str(jj)]);
    
    subplot(224);
    imagesc(angle(squeeze(Pmrho(jj,:,:))));
    colorbar;
    axis image;
    title(['Irene IFT( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    pause(.5); 

end

figure;

jj= 33;%49;%round(size(rhobest,3)/2);

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
 plot(thscantoplot,rock_Psi2,'*r'); 
 hold on; 
 plot(thscantoplot,rock_test,'ob');
 
 
 figure;
 
 for jj=1:size(Psi2,3) 
    
    
    subplot(221);
    imagesc(abs(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 fftn(rhobest) ' num2str(jj)]);
    
    subplot(222);imagesc(angle(Psi2_test(:,:,jj)));
    colorbar;
    axis image;
    title(['Psi2 fftn(rhobest) phase val ' num2str(jj)]);

    subplot(223);
    imagesc(abs(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) ) absolute val ' num2str(jj)]);
    
    subplot(224);
    imagesc(angle(Psi2(:,:,jj)./sqrt(size(Psi2,1)*size(Psi2,2))));
    colorbar;
    axis image;
    title(['Irene ( FT (rhobest) ) rhobest phase val ' num2str(jj)]); 
    
    pause(.5); 

 end


 
 figure;
 subplot(121);
 di(Psi2_test,-0.01);
 axis image
 title('FT conventional');
 
 subplot(122);
 di(Psi2,-0.01);
 axis image;
 title('FT Irene');