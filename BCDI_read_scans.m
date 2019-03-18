function [data_exp,rock,imgs] = BCDI_read_scans(filename,directname,specscan)

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

    % create a data_exp structure out of the tif files and the 
    figure;
    
    for ii=0:specscan.size-1 
         im = double(imread([directname filename num2str(ii,'%5.5i') '.tif'])); 
         imgs(:,:,ii+1)=im;
         
         subplot(121);
         imagesc(imgs(:,:,ii+1));
         axis image;
         title(num2str(ii+1));
         pause(.5);
         
         
         data_exp(ii+1).dth = thscan_list(ii+1);
         
         Ry = [cosd(-data_exp(ii+1).dth) 0 sind(-data_exp(ii+1).dth);
            0 1 0;
            -sind(-data_exp(ii+1).dth) 0 cosd(- data_exp(ii+1).dth)];
        
        ki = (Ry * ki_o.').';
        kf = (Ry * kf_o.').';
        
        subplot(122);
        quiver3(0,0,0, ki_o(1), ki_o(2), ki_o(3), 'r');
        hold on;
        quiver3(0,0,0, ki_o(1), ki_o(2), ki_o(3), 'r');
        quiver3(0,0,0, ki(1), ki(2), ki(3), 'k');
        quiver3(0,0,0, kf(1), kf(2), kf(3), 'k');
        
        dq_shift_nominal(ii+1,:) = (kf-ki)-qbragg;
        
        data_exp(ii+1).dqx = dq_shift_nominal(ii+1,1);
        data_exp(ii+1).dqy = dq_shift_nominal(ii+1,2);
        data_exp(ii+1).dqz = dq_shift_nominal(ii+1,3);
    
    
         
         data_exp(ii+1).I = imgs(:,:,ii+1);

    end

  
    
    % calculate the rocking curve:
    rock = squeeze(sum(sum(imgs),2));
    
    

end

