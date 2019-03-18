% this scripts replaces the measured intensities by calculated intensities
% in the case we want to use the simulated sample for the phase retrieval:

%{
%if(randomIniGuess)
    for ii = 1:numel(data_exp)
        mn = mean(data_exp(1).simI(:));
        data_exp(ii).simI = data_exp(ii).simI./mn;
    end
%end
%}

if(usesimI)
    display(['overwriting experimental I with sim I']);
    for ii=1:numel(data_exp)
        data_exp(ii).I = data_exp(ii).simI;
    end
end



figure(7); 
clf; setfigsize(gcf, 800,400);
colormap jetvar;
hax1 = axes('position', [0 0 .5 1]);
hax2 = axes('position', [0.5 0 .5 1]);
ca = [0 2];

for ii=1:numel(data_exp)
    %imagesc(hax1, log10(data_exp(ii).I));
    imagesc(hax1, (data_exp(ii).I));
    %caxis(hax1, ca); 
    set(hax1, 'xtick', Npix/2+1, 'ytick', Npix/2+1); %grid(hax1, 'on');
    drawnow;
    %imagesc(hax2, log10(data_exp(ii).simI)); axis image;
    imagesc(hax2, (data_exp(ii).simI)); axis image;
    %caxis(hax2, ca); 
    set(hax2, 'xtick', Npix/2+1, 'ytick', Npix/2+1); %grid(hax2, 'on');
    %caxis(ca);
    display(['dth: ' num2str(data_exp(ii).dth) 'dth_delta:' num2str(data_exp(ii).dth_delta)]);
    pause(.1);
end
