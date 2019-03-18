function [data_exp,rock, imgs_center] = BCDI_read_center_pad_scans(filename,directname,specscan,mindata)
    % This function read the data in filename, and stores the experimentl
    % data in a structure: data_exp, but it also centers the rocking curve
    % (the integrated intensity of each frame) by associating to the
    % maximum intensity of the rocking curve (the frame corresponding to
    % the Bragg peak) to the zero delta_theta. it also zero pads the
    % experimental data.
    
    global   ki_o kf_o
    
    % position of the bragg peak:
    qbragg = kf_o -ki_o;
    
    % read the experimental theta angles out of the specscan:
    %key= 'th';
    %start_index = findstr(specscan.scanline,'th');
    %str = specscan.scanline(start_index+strlength(key):end);
    %thscan_var = sscanf(str,'%f%f%d%d');
    %thscan_range = thscan_var(2) - thscan_var(1);
    thscan_list = specscan.var1;%[thscan_var(1):thscan_range/thscan_var(3):thscan_var(2)];
    %thscan_diff = thscan_list(2)-thscan_list(1);
    % create a data_exp structure out of the tif files and the 
    figure;
    
    for ii=0:specscan.size-1 
         im = double(imread([directname filename num2str(ii,'%5.5i') '.tif'])); 
         imgs(:,:,ii+1)=im;
    end
    
    [imgs_center, xyz,max_index] = center_array_Irene(imgs);
    
     % padding the data matrix to 
    arrysize = size(imgs_center,1);
    arrysize3=128;%size(imgs_center,3);%

    datap = zeros(arrysize,arrysize,arrysize3);

    datap(arrysize/2-size(imgs_center,1)/2+1:arrysize/2+size(imgs_center,1)/2,arrysize/2-size(imgs_center,1)/2+1:arrysize/2+size(imgs_center,1)/2, ...
    arrysize3/2-size(imgs_center,3)/2+1:arrysize3/2+size(imgs_center,3)/2) = imgs_center; %sample to the losless array size

    imgs_center = datap;
    
    ind=( imgs_center < mindata );
    imgs_center(ind)=0;
    
    % calculate the rocking curve:
    rock = squeeze(sum(sum(imgs_center),2));
    
    [val,index] = max(rock);
    
    thscan_padded(round(numel(rock)/2-specscan.size/2)+1:round(numel(rock)/2+specscan.size/2)) = thscan_list;

    
    thscan_list_center = thscan_padded-thscan_list(index);
    
    thscan_padded = 100.*ones(numel(rock),1);
    
    
%     thscan_diff = diff(thscan_list_center);
%     
%     thscan_diff_padded = zeros(numel(rock)-1,1);
%     
%     thscan_dif_padded(round((numel(rock)-1)/2-specscan.size/2)+1:round((numel(rock)-1)/2+specscan.size/2)) = thscan_diff;
%     
    for jj = find(thscan_padded == 100)' 
        if jj<round(size(thscan_padded,1)/2)
            first_nonzero_index = round(numel(rock)/2-specscan.size/2)+1;
            thscan_padded(jj) =  (jj-first_nonzero_index)*(thscan_list_center(first_nonzero_index+1)-thscan_list_center(first_nonzero_index))+thscan_list_center(1); 
        else
            last_nonzero_index = round(numel(rock)/2+specscan.size/2);
            thscan_padded(jj) =  (jj-last_nonzero_index)*(thscan_list_center(end)-thscan_list_center(end-1))+thscan_list_center(end); 
        end
    end
     
    thscan_list_center = thscan_padded;
    
    for ii=1:numel(thscan_list_center)
%          subplot(121);
%          imagesc(imgs(:,:,ii));
%          axis image;
%          title(['non centered sequence' num2str(ii)]);
         
         %subplot(122);
         imagesc(imgs_center(:,:,ii));
         axis image;
         title(['centered sequence' num2str(ii)]);
         pause(.5);
         
                 
         data_exp(ii).dth = thscan_list_center(ii);
         
         Ry = [cosd(-data_exp(ii).dth) 0 sind(-data_exp(ii).dth);
            0 1 0;
            -sind(-data_exp(ii).dth) 0 cosd(- data_exp(ii).dth)];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
        
%         subplot(122);
%         quiver3(0,0,0, ki_o(1), ki_o(2), ki_o(3), 'r');
%         hold on;
%         quiver3(0,0,0, ki_o(1), ki_o(2), ki_o(3), 'r');
%         quiver3(0,0,0, ki(1), ki(2), ki(3), 'k');
%         quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
%         
        dq_shift_nominal(ii,:) = (kf-ki)-qbragg;
        
        data_exp(ii).dqx = dq_shift_nominal(ii,1);
        data_exp(ii).dqy = dq_shift_nominal(ii,2);
        data_exp(ii).dqz = dq_shift_nominal(ii,3);
    
    
         
         data_exp(ii).I = sqrt(imgs_center(:,:,ii));

    end

  
end

