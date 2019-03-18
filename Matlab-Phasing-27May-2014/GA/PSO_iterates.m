function [pn_g velocity] = PSO_iterates(pn_g,pn_best,velocity,group_best,chi,params)

%breeds iterates based on selection criteria
%assume chi is array of chi values
%pn_g is the population of iterates

try
    psi2=params.psi2;
catch
    psi2=1.0;
end

dims=ndims(pn_g);
npop=size(pn_g,dims);

alpha=zero_phase(group_best);


amp_only=0;

pso_type='PSO';
%pso_type='FIPSO';

%
if amp_only == 1,disp('------AMP ONLY------'),end

phi1=2.05;
phi2=2.05;%

psi1=0.7298*psi2;

switch pso_type
    
    case 'PSO'
        
        for qq=1:npop

            %U1=random('Uniform',0,1);  %1.5;%
            %U2=random('Uniform',0,1);  %1.5;%
            U1=random('Uniform',0,1,size(alpha));  %1.5;%
            U2=random('Uniform',0,1,size(alpha));  %1.5;%


        if dims == 3
            % align current posn with group best
            pn=align_and_zero(alpha,pn_g(:,:,qq));
            % align best current posn with group best
            pn_b=align_and_zero(alpha,pn_best(:,:,qq));
            Vold=velocity(:,:,qq);
        end
        if dims == 4
            % align current posn with group best
            pn=align_and_zero(alpha,pn_g(:,:,:,qq));
            % align best current posn with group best
            pn_b=align_and_zero(alpha,pn_best(:,:,:,qq));
            Vold=velocity(:,:,:,qq);
        end

        if amp_only == 1
            % personeal influence
            Pinf=abs(pn_b)-abs(pn);
            % group influence
            Ginf=abs(alpha)-abs(pn);
            Vnew=Vold+phi1.*U1.*Pinf+phi2.*U2.*Ginf;
        else
            % personeal influence
            Pinf=pn_b-pn;
            % group influence
            Ginf=alpha-pn;
            Vnew=Vold+phi1.*U1.*Pinf+phi2.*U2.*Ginf;
        end

        % move to the new position
        if dims == 3,
            velocity(:,:,qq)=Vnew;
            pn_g(:,:,qq)=pn+psi1*Vnew;
        end

        if dims == 4,
            velocity(:,:,:,qq)=Vnew;
            pn_g(:,:,:,qq)=pn+psi1*Vnew;
        end


        end

    case 'FIPSO'
        
        %get the order of bestness
        [chi_ordered chi_index_l]=sort(chi);
        n_neighbours=3;
        chi_index=chi_index_l(1:n_neighbours);
        
        disp(['------USING - [',num2str(numel(chi_index)),'] NEIGHBOURS------'])
        
        for qq=1:npop

            if dims == 3
                pn=zero_phase(pn_g(:,:,qq)); %zero the current posn
                Vold=velocity(:,:,qq);
            end

            if dims == 4
                pn=zero_phase(pn_g(:,:,:,qq));
                Vold=velocity(:,:,:,qq);
            end
                    
            pn_x=0;                     %get the array for the sum of differences
            
            %loop through all to get the differences relative to the current
            %only include contribution from the X best, given by chi_index

            for ww=1:numel(chi_index)               
                
                bb=chi_index(ww);
                
                if dims == 3    %zero the current best and align with the current
                    pn_b=align_and_zero(pn,pn_best(:,:,bb));
                end
                if dims == 4
                    pn_b=align_and_zero(pn,pn_best(:,:,:,bb));
                end

                U1=random('Uniform',0,1,[size(alpha)]);  %random arrays

                pn_x=pn_x+U1.*(pn_b-pn);     %the sum of differences between current posn and all the best

            end
            
            pn_x=pn_x/npop; %normalize with repsect to the npop
            
            Vnew=Vold+phi1.*pn_x;

            % move to the new position
            if dims == 3,
                velocity(:,:,qq)=Vnew;
                pn_g(:,:,qq)=pn+psi1*Vnew;
            end

            if dims == 4,
                velocity(:,:,:,qq)=Vnew;
                pn_g(:,:,:,qq)=pn+psi1*Vnew;
            end
            
            
            
        end

        
end


end

function pn=align_and_zero(ref,pn)

%check the conj reflection
pn=check_conj_ref(ref,pn);

%align just the amplitudes intiailly
%shift iterates to the dominant one
[h k l]=register_3d_reconstruction(abs(ref),abs(pn));
pn=zero_phase(((sub_pixel_shift((pn),h,k,l))));
%pn=zero_phase(pn);


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

function pn=remove_ramp(pn)

amp=abs(pn);
ph=atan2(imag(pn),real(pn));

F0=ifftshift(fftn(fftshift(amp)));
F1=ifftshift(fftn(fftshift(amp.*exp(i*ph))));
[h k l]=register_3d_reconstruction(abs(F0),abs(F1));
pn=ifftshift(ifftn(fftshift(sub_pixel_shift(F1,h,k,l))));

end

function pn=check_conj_ref(ref,pn)

cnj_rf=is_conj_ref(abs(ref),abs(pn));
            
if cnj_rf ~=0,pn=conj_reflect(pn);end

end

function cnj_rf=is_conj_ref(a,b)

c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));

%c1=cross_correlation(a,b);
%c2=cross_correlation(a,b);


c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end

function [shifted]=sub_pixel_shift(array,row_shift,col_shift,z_shift)

    buf2ft=fftn(array);
    [nr,nc,nz]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    Nz = ifftshift([-fix(nz/2):ceil(nz/2)-1]);
    [Nc,Nr,Nz] = meshgrid(Nc,Nr,Nz);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc-z_shift*Nz/nz));
    %Greg = Greg*exp(i*diffphase);
    shifted=ifftn(Greg);

end