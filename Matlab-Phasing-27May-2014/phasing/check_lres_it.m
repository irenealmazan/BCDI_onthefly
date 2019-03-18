function [ params ] = check_lres_it(params,flip,qq)
%jclark
%check if the lres should be enabled, just set the gauss width to 1 so it
%will skip it

%params.sigma=params.sw_sig_orig;

%ALG1
if flip == 1

    if params.nonGA_lres_ALG1 == 0
        params.nonGA_lres_det(qq)=1;
        params.nonGA_lres_sig(qq)=params.sw_sig_orig;    
    end
    
else  %ALG2
    
    if params.nonGA_lres_ALG2 == 0
        params.nonGA_lres_det(qq)=1;
        params.nonGA_lres_sig(qq)=params.sw_sig_orig;
    end
    
end

params.sigma=params.nonGA_lres_sig(qq);

end

