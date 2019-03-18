function [mod_ini_square,mod_fin_square] = NW_calculate_final_dq_improvement(data_exp)

% this function calculates the distance between the initial dq_shift, the
% final dq_shift and the true dq_shift

    dq_ini = zeros(numel(data_exp),3);
    mod_ini_square = zeros(numel(data_exp),1);
    mod_fin_square = zeros(numel(data_exp),1);

    for ii = 1:numel(data_exp)  
        dq_ini(ii,:) = [data_exp(ii).dqx data_exp(ii).dqy data_exp(ii).dqz] ;
        mod_ini_square(ii) = (data_exp(ii).dqshift_delta(1) - dq_ini(ii,1))^2 +  (data_exp(ii).dqshift_delta(2) - dq_ini(ii,2))^2 + (data_exp(ii).dqshift_delta(3) - dq_ini(ii,3))^2;
        mod_fin_square(ii) = (data_exp(ii).dqshift_delta(1) - data_exp(ii).dqshift(1))^2 +  (data_exp(ii).dqshift_delta(2) - data_exp(ii).dqshift(2))^2 + (data_exp(ii).dqshift_delta(3) - data_exp(ii).dqshift(3))^2;
    end

end

