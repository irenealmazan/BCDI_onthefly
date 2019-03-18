function [ params ] = check_ALG_ptych( params )
%jclark
if params.flip == 1
    params.ALG=params.ALG1;
else
    params.ALG=params.ALG2;
end


end

