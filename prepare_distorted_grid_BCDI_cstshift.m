function [dq_shift,dth_delta_list] = prepare_distorted_grid_BCDI_cstshift(dth_disp,fly2Danglist,thBragg,index_to_distort,Ry2,Rx)

    global ki_o kf_o
    
  
    qbragg = kf_o - ki_o;
    
    
    dq_shift = zeros(numel(fly2Danglist),3);
    shift_delta = zeros(numel(fly2Danglist),3);

    
    dth_old = 1000;
    
    %index = randi(numel(dth_disp));
    
    
    for ii = 1:numel(fly2Danglist)
        
        dth_nominal = fly2Danglist(ii)-thBragg;
       % for jj = 1:numel(index_to_distort)
            if ~isempty(find(ii == index_to_distort))
                
                index = ii;
                
                dth_delta = dth_disp(index);
                
                dth = dth_nominal + dth_delta;
                
                Ry = [cosd(-dth) 0 sind(-dth);
                    0 1 0;
                    -sind(-dth) 0 cosd(-dth)];
                
                ki = (Ry * ki_o.').';
                kf = (Ry * kf_o.').';
                dq_shift(ii,:) = (kf-ki)-qbragg;
                dth_delta_list(ii) = dth_delta;
                
                
            else
                dth_delta = 0;
                dth = dth_nominal + dth_delta;
                
                Ry = [cosd(-dth) 0 sind(-dth);
                    0 1 0;
                    -sind(-dth) 0 cosd(-dth)];
                
                ki = (Ry * ki_o.').';
                kf = (Ry * kf_o.').';
                dq_shift(ii,:) = (kf-ki)-qbragg;
                dth_delta_list(ii) = dth_delta;
                
            end
     %   end
        
        
        
        
    end
    
  
    
    %{
    switch flag

        case 'theta'
            %%% DISTORT THETA GRID
            %{
           
            % add a random constant
            %index = randi(numel(dth_disp));
            %dth_delta = dth_disp(index);
            dth_old = 1000;
            
            for ii = 1:numel(data_exp)                               
                
                dth_nominal = data_exp(ii).dth;
                
                if dth_nominal == dth_old
                    index = randi(numel(dth_disp));
                    dth_delta = dth_disp(index);
                end
                
                dth = dth_nominal + dth_delta;
                
                Ry = [cosd(-dth) 0 sind(-dth);
                    0 1 0;
                    -sind(-dth) 0 cosd(-dth)];
                
                ki = (Ry * ki_o.').';
                kf = (Ry * kf_o.').';
                dq_shift(ii,:) = (kf-ki)-qbragg;
                
                dth_old = dth_nominal;
            end


            %}

            shift_delta = [0 0 0];

        case 'xy'
            
            for ii = 1:numel(data_exp)
                
                dth_nominal = data_exp(ii).dth;
                
                if dth_nominal == dth_old
                    index = randi(numel(dth_disp));
                    dth_delta = dth_disp(index);
                end
                
                dth = dth_nominal + dth_delta;
                
                Ry = [cosd(-dth) 0 sind(-dth);
                    0 1 0;
                    -sind(-dth) 0 cosd(-dth)];
                
                ki = (Ry * ki_o.').';
                kf = (Ry * kf_o.').';
                dq_shift(ii,:) = (kf-ki)-qbragg;
                
                dth_old = dth_nominal;
                
                % rotation matrices:
                
                Ry1=[cosd(th) 0 sind(th);
                    0 1 0;
                    -sind(th) 0 cosd(th)];
                
                %%%% DISTORT POSITION GRID
                %{
                index = randi(numel(dx_disp));
                xshift_delta = dx_disp(index);
                index = randi(numel(dy_disp));
                yshift_delta = dy_disp(index);
                index = randi(numel(dz_disp));
                zshift_delta = dz_disp(index); % note that a shift along one of the directions of the scan plane, involves a shift in x and z
                
                % build the shift vector, to lie in the scan plane:
                shift_delta(ii,:) = [xshift_delta yshift_delta zshift_delta];
                
                shift_delta(ii,:) = (Ry1*shift_delta')';
                shift_delta(ii,:) = (Ry2*shift_delta')';
                shift_delta(ii,:) = (Rx*shift_delta')';
                
                %}
                
                dq_shift(ii,:) = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz];
                
            end
        case 'both'
            dth_nominal = data_exp.dth;

            % add a random constant
            index = randi(numel(dth_disp));
            dth_delta = dth_disp(index);

            dth = dth_nominal + dth_delta;

            Ry = [cosd(-dth) 0 sind(-dth);
                0 1 0;
                -sind(-dth) 0 cosd(-dth)];

            ki = (Ry * ki_o.').';
            kf = (Ry * kf_o.').';

            dq_shift = (kf-ki)-qbragg;

            %%%% DISTORT POSITION GRID
            %%{
            index = randi(numel(dx_disp));
            xshift_delta = dx_disp(index);
            index = randi(numel(dy_disp));
            yshift_delta = dy_disp(index);
            index = randi(numel(dz_disp));
            zshift_delta = dz_disp(index); % note that a shift along one of the directions of the scan plane, involves a shift in x and z

            % build the shift vector, to lie in the scan plane:
            shift_delta = [xshift_delta yshift_delta zshift_delta];

            shift_delta = (Ry1*shift_delta')';
            shift_delta = (Ry2*shift_delta')';
            shift_delta = (Rx*shift_delta')';



    end

%}
end