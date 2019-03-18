function [cx_its] = average_iterates( iterates)
%UNTITLED2 Summary of this function goes here
%   
% calcs an average from iterates

dims=ndims(iterates);

if dims == 3, alpha = (iterates(:,:,1));end
if dims == 4, alpha = (iterates(:,:,:,1));end

avg_its=abs(alpha);

max_its=abs(alpha);

cx_its=zero_phase(alpha);

n_its=size(iterates,dims);

for qq = 2:n_its

    if dims == 3, beta = (iterates(:,:,qq));end
    if dims == 4, beta = (iterates(:,:,:,qq));end
    
    beta=check_conj_ref(alpha,beta);

    %align just the amplitudes intiailly
    %shift the alpha TO the beta
    [h k l]=register_3d_reconstruction(abs(alpha),abs(beta));
    
    beta_s=((sub_pixel_shift((beta),h,k,l)));
    
    avg_its=avg_its+abs(beta_s);

    max_its=max(max_its,abs(beta_s));
    
    cx_its=cx_its+zero_phase(beta_s);
    
    beta=[];
    beta_s=[];
    
end

avg_its=avg_its/dims;


end

function pn=zero_phase(pn,val)

try
    val;
catch
    val=0;
end

ph=atan2(imag(pn),real(pn));

SS=shrink_wrap(abs(pn),.2,.5);  %get just the crystal, i.e very tight support
avg_ph=mean(ph(find(SS > 0)));
ph=ph-avg_ph+val;

pn=abs(pn).*exp(i*ph);


end