
clear data_sim data_exp

%%%% select the data_exp index which are of interest:
%{
dth_old = 1000;
cnt = 1;
for ii = 1:numel(data_exp)
    dth_nominal = data_exp(ii).dth;
        if dth_nominal ~= dth_old
            index(cnt) = ii;
            cnt = cnt + 1;
        end
   dth_old = dth_nominal;     
end

%%% rewrite data_exp for bcdi simulation
data_exp_new = data_exp;

clear data_exp
        
for ii = 1:numel(index)
   data_exp(ii) = data_exp_new(index(ii)); 
end
%}

global ki_o ki_f

 qbragg = kf_o - ki_o;

mxI = zeros(size(thscanlist));
im_sum_sim = zeros(Npix);

ki_norm = ki_o/norm(ki_o);
kf_norm = [0 0 1];


angx = acosd(dot(ki_norm([1 3]), kf_norm([1 3])));
angy = acosd(dot(ki_norm([2 3]), kf_norm([2 3])));



Ry2=[cosd(-del) 0 sind(-del);
    0 1 0;
    -sind(-del) 0 cosd(-del)];

Rx=[1 0 0;
    0 cosd(gam) -sind(gam);
    0 sind(gam) cosd(gam)];

%%%% NEW: generate random displacements of the beam with respect to its ideal
% position: 
dx_disp = [0]*d2_bragg;% [-10:5:10]*d2_bragg;
dy_disp = [0]*d2_bragg;
dz_disp = [0]*d2_bragg;
 %index = [1 2 2];
%

% theta:


index_to_distort = [34];%[28:1:45];%[28 31 35 36];%[36];%[28:1:33];%randi(numel(fly2Danglist),1,number_angles_distort);

dth_disp(index_to_distort) = [0.02];%[0.008];%[0.017 -0.008 -0.013 0.005];%[0.017];% -0.008 -0.013 0.005];%[0.002];%[-.003:.001:.003];

[dq_shift,dth_delta_list] = prepare_distorted_grid_BCDI_cstshift(dth_disp,fly2Danglist,thBragg,index_to_distort,Ry2,Rx);

rock_curve = zeros(numel(fly2Danglist),1);

figure(27);clf;
for ii = 1:numel(fly2Danglist)
    
    
    % build the nominal angle grid:
    data_exp(ii).dth = fly2Danglist(ii)-thBragg;% + dth_delta_list(ii);
    
    Ry = [cosd(-data_exp(ii).dth) 0 sind(-data_exp(ii).dth);
            0 1 0;
            -sind(-data_exp(ii).dth) 0 cosd(- data_exp(ii).dth)];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
        dq_shift_nominal(ii,:) = (kf-ki)-qbragg;
        
        data_exp(ii).dqx = dq_shift_nominal(ii,1);
        data_exp(ii).dqy = dq_shift_nominal(ii,2);
        data_exp(ii).dqz = dq_shift_nominal(ii,3);
    
    % calculate the Q term which depends on the theta angle distorted!!! :
     Qterm = exp(i* dq_shift(ii,1) * X) .* ...
            exp(i* dq_shift(ii,2) * Y) .* ...
            exp(i* dq_shift(ii,3) * Z);    
    
     
    
    
    data_exp(ii).dqshift_delta = dq_shift(ii,:);
    data_exp(ii).dth_delta = dth_delta_list(ii);
    
    tempbm = probe;
    %tempbm = circshift(probe, round([yshift_shifted  xshift_shifted zshift_shifted]/d2_bragg));
    
        %{
        clf; hold on;
        %h=di(temp,-.5,'y',X,Y,Z); alpha(h,.5); 
        h=di(NW, -.5, 'r', X,Y,Z); alpha(h, .2); 
        h=di(tempbm, -.5, 'c', X,Y,Z); alpha(h, .2);
        %plot3(pos_probe_exp(:,1), pos_probe_exp( :,2), pos_probe_exp(:,3), 'k.');
        %scatter3(pos_probe_exp(:,1), pos_probe_exp(:,2), pos_probe_exp(:,3), 50, Ga_chan, 'filled');
        plot3(data_exp(ii).xpos, data_exp(ii).ypos, data_exp(ii).zpos, 'ro');
        %view([28.4, 52.2]);
        view([-11 -6]);
        pause(.1); hold off;
        %return
        %}
    
    

    if(addNWsf)
        temp = sum( NWsfzb.*tempbm.*Qterm, 3);
    else
        temp = sum( NW.*tempbm.*Qterm, 3);
    end

    if(mod(ii,1)==0) 
        display(['simulating dp, ' num2str(ii) ' of ' num2str(numel(data_exp))]); 
        subplot(121); 
        imagecomp(temp); axis image; drawnow;
    end

    %data_exp(ii).simproj = temp;
    %data_exp(ii).simamp = fftshift(fftn(fftshift(data_exp(ii).simproj)));
    temp = fftshift(fftn(fftshift(temp)));
    temp = flipud(temp);
    
    if(mod(ii,1)==0) subplot(122); imagesc(abs(temp));axis image; drawnow; end
    
    data_exp(ii).simI = sqrt(temp .* conj(temp));
    mxI(ii) = max(data_exp(ii).simI(:));
    
    im_sum_sim = im_sum_sim + data_exp(ii).simI;
   
    % the 3D bragg peak
    Psi2_toplot(:,:,ii) = temp;
    
    rock_curve(ii) = sum(sum(data_exp(ii).simI,1));
    
end

middpind = round(numel(data_exp)/2);

    
figure;
plot(fly2Danglist-thBragg,rock_curve,'*r');
title('distorted rocking curve');