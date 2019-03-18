function [ output_args ] = batch_output_convert_rec(top_dir,prefix,prefix2)
%jclark
%do a batch conversion from rec frame to lab frame
if exist('prefix2') == 0
    prefix2='';
end

string=[prefix,'**-AMP.rec'];
amp_fs=rdir([top_dir,'**/*/',prefix2,string]);

string=[prefix,'**-PARAMS.mat'];
pm_fs=rdir([top_dir,'**/*/',prefix2,string]);

string=[prefix,'**-PH.rec'];
ph_fs=rdir([top_dir,'**/*/',prefix2,string]);

n_temps=size(amp_fs,1);


for qq=1:n_temps
    
    load(pm_fs(1).name)
    ddir=extract_dir_from_string(amp_fs(qq).name);
    output_convert_rec(ddir)

end



end

