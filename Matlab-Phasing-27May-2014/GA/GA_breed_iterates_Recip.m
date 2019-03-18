function [ pn_g ] = GA_breed_iterates_Recip(pn_g,chi,params,best,data)

%breeds iterates based on selection criteria
%assume chi is array of chi values
%pn_g is the population of iterates

% a is the 'best' iterate and is kept constant for each other iterate, b
% sqrt_ab -  traditional GHIO,sqrt(a*b)
% max_ab  -  my take, max(|a|,|b|)exp(1/2(pa+pb))
% avg_ab  -  my take, 1/2a+1/2b

% sqrt_ab_pa - sqrt(|a||b|)exp(pa)
% max_ab_pa  - max(|a|,|b|)exp(pa)
% avg_ab_pa  - 1/2(|a|+|b|)exp(pa)

% sqrt_abg -  (a*b*g)^(1/3)
% max_abg  -  max(|a|,|b|,|g|)exp(1/3(pa+pb+pb))
% avg_abg  -  1/3(a+b+g)

% sqrt_abg_pa - (|a||b||g|)^(1/3)exp(pa)
% max_abg_pa  - max(|a|,|b|,|g|)exp(pa)
% avg_abg_pa  - 1/2(|a|+|b|+|g|)exp(pa)


val=0;
dims=ndims(pn_g);
npop=size(pn_g,dims);

%get the worst one, will replace it
ind_b=find(chi ==max(chi(:)));

%get the 'alpha male'
ind=find(chi == min(chi(:)));

if numel(ind) > 1,ind=ind(1);end

if dims == 3, alpha=pn_g(:,:,ind);end
if dims == 4, alpha=pn_g(:,:,:,ind);end

%calc the average for use in acg_abc
% if strcmp(params.breed_mode,'avg_abg') == 1,
%     aligned=abs(align_iterates(pn_g,ind));
%     avg_amp=mean(aligned,dims);           %align with the alpha
%     sig_amp=std(aligned,1,dims);
% end

try
    breed_mode=params.breed_mode;
catch
    breed_mode='sqrt_ab';
end


%get the gamma male, the second best
chi_g=chi;%
chi_g(ind)=max(chi(:)); %do this to get the second one
ind_g=find(chi_g == min(chi_g(:)));

if numel(ind_g) > 1,ind_g=ind_g(1);end

if dims == 3, gamma=pn_g(:,:,ind_g);end
if dims == 4, gamma=pn_g(:,:,:,ind_g);end
gamma=zero_phase(gamma,val);
ph_gamma=atan2(imag(gamma),real(gamma));


alpha_avg=0;
alpha=zero_phase(alpha,val);            

x0=[1,1,1]*(1/3);  %fir use in the min sharp


for qq=1:npop
    
    if dims == 3, beta=pn_g(:,:,qq);end
    if dims == 4, beta=pn_g(:,:,:,qq);end
    
    if qq == 1,max_amp=max(abs(pn_g),[],3);end
    
    if qq ~= ind
        
        beta=zero_phase(beta,val);             
        
        %%ph_beta=atan2(imag(beta),real(beta));

        %align just the amplitudes intiailly
        %shift the alpha TO the beta
        %check the conj reflection
        
        alpha=check_conj_ref(beta,alpha);
        %alpha_s=zero_phase(alpha,val);
        
        [h k l]=register_3d_reconstruction(abs(beta),abs(alpha));
        alpha_s=zero_phase((sub_pixel_shift(alpha,h,k,l)),val);
        %%ph_alpha=atan2(imag(alpha_s),real(alpha_s));
        
        gamma=check_conj_ref(beta,gamma);
        %gamma_s=zero_phase(gamma,val);
        
        [h k l]=register_3d_reconstruction(abs(beta),abs(gamma));
        gamma_s=zero_phase((sub_pixel_shift(gamma,h,k,l)),val);
        %%ph_gamma=atan2(imag(gamma_s),real(gamma_s));
        
        Fbeta=fftshift(fftn(fftshift(beta)));
        ph_beta=atan2(imag(Fbeta),real(Fbeta));
        
        Fgamma=fftshift(fftn(fftshift(gamma_s)));
        ph_gamma=atan2(imag(Fgamma),real(Fgamma));
        
        Falpha=fftshift(fftn(fftshift(alpha)));
        ph_alpha=atan2(imag(Falpha),real(Falpha));
        
        switch breed_mode
                
            case {'sqrt_ab'}
                beta_t=(abs(Fbeta).*abs(Fgamma)).^.5.*exp(i*0.5*(ph_beta+ph_gamma));
                beta_t=fftshift(ifftn(fftshift(beta_t))); 
            case {'max_ab'}
                
                
            
                
        end
        

        if dims == 3, pn_g(:,:,qq)=beta_t;end
        if dims == 4, pn_g(:,:,:,qq)=beta_t;end

        beta=[];
        beta_t=[];
        ph_beta=[];
        alpha_s=[];
 
    end
        
end



beta=[];
beta_t=[];
ph_beta=[];
alpha_s=[];
max_amp=[];

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

%c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
%c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));

c1=cross_correlation(a,b);
c2=cross_correlation(a,b);


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

