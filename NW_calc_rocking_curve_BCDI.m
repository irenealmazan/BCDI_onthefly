
% clear rock_probe rock_pb dqlist_rock
% 
% thscanlist_rock = fly2Danglist-thBragg+dth_delta_list;%[-.8:.01:.8];
% 
% dqlist_rock = zeros(numel(thscanlist_rock), 3);

for ii=1:numel(data_exp)
    
    thscan = data_exp(ii).dth;

    Ry = [cosd(thscan) 0 sind(thscan);
        0 1 0;
        -sind(thscan) 0 cosd(thscan)];
    ki = (Ry *ki_o')';
    kf = (Ry *kf_o')';
    
    dqlist_rock(ii,:) = (kf-ki) - qbragg;

end    
    
%%

cenx_rock = [];
ceny_rock = [];
ax_temp = [-Npix/2:Npix/2-1];
ax_temp = ax_temp(:);

for ii=1:numel(thscanlist_rock)
    
    Qterm = exp(i* dqlist_rock(ii,1) * X) .* ...
            exp(i* dqlist_rock(ii,2) * Y) .* ...
            exp(i* dqlist_rock(ii,3) * Z);
        
    temp_probe = sum( (NW).*probe.*Qterm, 3);
    temp_pb = sum( (NW).*Qterm, 3);
    
    temp_dp_probe = fftshift(fftn(fftshift(temp_probe)));
    
    rock_probe(:,:,ii) = temp_dp_probe;
    rock_pb(:,:,ii) = fftshift(fftn(fftshift(temp_pb)));        
   
    temp = sum(temp_dp_probe.*conj(temp_dp_probe), 1);
    temp = temp(:);
    cenx_rock = [cenx_rock sum(ax_temp.*temp)/sum(temp)];

    temp = sum(temp_dp_probe.*conj(temp_dp_probe), 2);
    temp = temp(:);
    ceny_rock = [ceny_rock sum(ax_temp.*temp)/sum(temp)];

end



figure(4);
clf; 
subplot(221); di(rock_pb,-.05);
subplot(222); plot(thscanlist_rock, squeeze(sum(sum(abs(rock_pb),2),1)));
xlabel('rock angle'); title('unfocused beam'); set(gca, 'xtick', thscanlist, 'XGrid', 'on');
subplot(223); di(rock_probe,-.05);
subplot(224); plot(thscanlist_rock, squeeze(sum(sum(abs(rock_probe),2),1)));
xlabel('rock angle'); title('nano-focused beam'); set(gca, 'xtick', thscanlist, 'xgrid','on');
drawnow

