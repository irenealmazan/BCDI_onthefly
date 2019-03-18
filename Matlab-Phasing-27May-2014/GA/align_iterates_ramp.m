function [aligned_its] = align_iterates_ramp( iterates,ind,val)
%Jclark  
% aligns a whole series of iterates
% ind is the one to align to. default == 1

try
    ind;
catch
    ind = 1;
end
try
    val;
catch
    val=[];
end

dims=ndims(iterates);

if dims == 3, alpha = (iterates(:,:,ind));end
if dims == 4, alpha = (iterates(:,:,:,ind));end

alpha=remove_ramp_pn_ups(alpha,3);

aligned_its=zeros(size(iterates));
n_its=size(iterates,dims);


for qq = 1:n_its

    
    if dims == 3, beta = (iterates(:,:,qq));end
    if dims == 4, beta = (iterates(:,:,:,qq));end
    
    if qq ~= ind
        beta=check_conj_ref(alpha,beta);

        beta=remove_ramp_pn_ups(beta,3);
        
        %align just the amplitudes intiailly
        %shift the alpha TO the beta
        [h k l]=register_3d_reconstruction(abs(alpha),abs(beta));

        beta_s=((sub_pixel_shift((beta),h,k,l)));
    else
        beta_s=beta;
    end
        
    if isempty(val) == 1
        if dims == 3, aligned_its(:,:,qq)=(beta_s);end
        if dims == 4, aligned_its(:,:,:,qq)=(beta_s);end
    else
        if dims == 3, aligned_its(:,:,qq)=zero_phase(beta_s,val);end
        if dims == 4, aligned_its(:,:,:,qq)=zero_phase(beta_s,val);end
    end
    
    beta=[];
    beta_s=[];
    
end


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