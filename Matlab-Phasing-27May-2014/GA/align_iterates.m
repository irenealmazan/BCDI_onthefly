function [aligned_its] = align_iterates( iterates,ind,val)
%Jclark  
% aligns a whole series of iterates
% ind is the one to align to. default == 1
% use ind = -1 to align sequenialy, ie. 2 -1,3-2,4-3 etc

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


aligned_its=zeros(size(iterates));
n_its=size(iterates,dims);

switch ind
    
    case -1

        for qq = 1:n_its

            
            if dims == 3, beta = (iterates(:,:,qq));end
            if dims == 4, beta = (iterates(:,:,:,qq));end

            if qq > 1
                if dims == 3, alpha = (iterates(:,:,qq-1));end
                if dims == 4, alpha = (iterates(:,:,:,qq-1));end
            else
                alpha=beta;
            end
            
            beta=check_conj_ref(alpha,beta);

            %align just the amplitudes intiailly
            %shift the alpha TO the beta
            [h k l]=register_3d_reconstruction(abs(alpha),abs(beta));

            beta_s=((sub_pixel_shift((beta),h,k,l)));

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
        
    otherwise
        
        if dims == 3, alpha = (iterates(:,:,ind));end
        if dims == 4, alpha = (iterates(:,:,:,ind));end
        
        for qq = 1:n_its


            if dims == 3, beta = (iterates(:,:,qq));end
            if dims == 4, beta = (iterates(:,:,:,qq));end

            if qq ~= ind
                beta=check_conj_ref(alpha,beta);

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