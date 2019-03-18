function [ params ] = check_lres(params)
%jclark

params.sw_sig_orig=params.sigma;  %store original

%check if alg1 has lres enabled
if strcmp('lr',params.ALG1(end-1:end)) 
    params.nonGA_lres_ALG1=1;
    disp(' ')
    disp('Reconstructing from low to high resolution for ALG1....')
    disp(' ')
else
    params.nonGA_lres_ALG1=0;
    params.nonGA_lres_sig=make_vals(params.iterations,params.sw_sig_orig,params.sw_sig_orig);
    params.nonGA_lres_det=make_vals(params.iterations,1,1);
end

%check if alg2 has lres enabled
if strcmp('lr',params.ALG2(end-1:end)) 
    params.nonGA_lres_ALG2=1;
    disp(' ')
    disp('Reconstructing from low to high resolution for ALG2....')
    disp(' ')
else
    params.nonGA_lres_ALG2=0;
    params.nonGA_lres_sig=make_vals(params.iterations,params.sw_sig_orig,params.sw_sig_orig);
    params.nonGA_lres_det=make_vals(params.iterations,1,1);
end

%sw sigma max
if isfield(params,'nonGA_sig_max') == 0
    params.nonGA_sig_max=3;
end

%with of gauss at the end
if isfield(params,'nonGA_lres_detfinal') == 0
    params.nonGA_lres_detfinal=1;
end

%width of gauss initially
if isfield(params,'nonGA_lres_detinit') == 0
    params.nonGA_lres_detinit=.1;
end

%sw values
if isfield(params,'nonGA_lres_sig') == 0
    params.nonGA_lres_sig=make_vals(params.iterations,params.sw_sig_orig,params.nonGA_sig_max);
else
    if numel(params.nonGA_lres_sig) ~= params.iterations,params.nonGA_lres_sig(end:params.iterations)=params.nonGA_lres_sig(end);end 
end

%gauss width values
if isfield(params,'nonGA_lres_det') == 0
    params.nonGA_lres_det=make_vals(params.iterations,params.nonGA_lres_detfinal,params.nonGA_lres_detinit);
else
    if numel(params.nonGA_lres_det) ~= params.iterations,params.nonGA_lres_det(end:params.iterations)=params.nonGA_lres_det(end);end 
end

end

