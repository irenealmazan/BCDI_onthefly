function [ output_args ] = batch_change_params_mat(top_dir,prefix,ex_string)
%jclark
%change  value in a params file. useful if you get the geometry wrong


string=[prefix,'**-PARAMS.mat'];
pm_fs=rdir([top_dir,'**/*/',string]);


n_temps=size(pm_fs,1);

for qq=1:n_temps
    
    load(pm_fs(qq).name)
    disp(ex_string)
    eval(ex_string)
    save(pm_fs(qq).name,'params')

end


end

