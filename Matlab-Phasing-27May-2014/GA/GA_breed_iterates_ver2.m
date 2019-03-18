function [ pn_g ] = GA_breed_iterates_ver2(pn_g,chi,params,best,data)

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

switch params.breed_mode

    case {'Dhalf-best'}
        disp('Using best (from all previous generations)....')
        if dims == 3, pn_g(:,:,ind)=best;end
        if dims == 4, pn_g(:,:,:,ind)=best;end
end

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

%get the order of best
[B IX]=sort(chi);

% align and average all iterates to get the mean and sig
switch breed_mode

    case {'N-dist','Dhalf','avg_all','max_all','orth'}
        [summed sigma pn_g]=sum_align_GA_reconstructions(pn_g,ind);
end

if strcmp(breed_mode,'orth') == 1
 
    [ pn_g ] = orth_iterates(pn_g);
    
else
    
        
    for qq=1:npop

        if dims == 3, beta=pn_g(:,:,qq);end
        if dims == 4, beta=pn_g(:,:,:,qq);end

        if qq == 1,max_amp=max(abs(pn_g),[],3);end

        if qq ~= ind

            beta=zero_phase(beta,val);             
            ph_beta=atan2(imag(beta),real(beta));

            %align just the amplitudes intiailly
            %shift the alpha TO the beta
            %check the conj reflection

            if strcmp(breed_mode,'N-dist') == 0   %no need since they are all aligned
                alpha=check_conj_ref(beta,alpha);
                [h k l]=register_3d_reconstruction(abs(beta),abs(alpha));

                alpha_s=zero_phase((sub_pixel_shift(alpha,h,k,l)),val);
                ph_alpha=atan2(imag(alpha_s),real(alpha_s));

                gamma=check_conj_ref(beta,gamma);
                [h k l]=register_3d_reconstruction(abs(beta),abs(gamma));
                gamma_s=zero_phase((sub_pixel_shift(gamma,h,k,l)),val);
                ph_gamma=atan2(imag(gamma_s),real(gamma_s));
            end

            switch breed_mode


                case {'auto_focus'}
                    beta_t=auto_focus(beta,10e-9,1.39e-10,10e-6,'x');

                case {'pixel_switch'}
                    beta_t = random_mix_array( beta,alpha_s,0.5 );

                case {'max_all'}
                    amp=max(abs(pn_g),[],dims);
                    beta_t=amp.*exp(i*ph_beta);

                case {'avg_all'}
                    beta_t=summed.*exp(i*ph_beta);

                case {'Dhalf','Dhalf-best'}
                    Nhalf=round(numel(IX)/2);
                    if dims == 3
                        Delta=Nhalf*pn_g(:,:,qq)-sum(pn_g(:,:,IX(1:Nhalf)),3);
                    end
                    if dims == 4
                        Delta=Nhalf*pn_g(:,:,:,qq)-sum(pn_g(:,:,:,IX(1:Nhalf)),4);
                    end
                    beta_t=beta+Delta;

                case {'dsqrt'}    
                    amp=abs(beta).^.5;
                    beta_t=amp.*exp(i*ph_beta);
                case {'N-dist'}
                    beta_t=beta+sigma;

                case {'avg_sqrt'}
                    bb=1/3;
                    amp=( (abs(beta)).^bb+(abs(alpha_s)).^bb+(abs(gamma_s)).^bb)/3;

                    amp=amp.^(1/bb);
                    beta_t=amp.* exp(i*(ph_beta));

                case {'b_pa'}
                    beta_t=abs(beta).* exp(i*(ph_alpha));

                case {'2ab_a_b'}
                    beta_t=2*(beta.*alpha_s)./(beta+alpha_s);

                case {'2a-b_pa'}
                    beta_t=(2*abs(alpha_s)-abs(beta)).* exp(i*(ph_alpha));

                case {'sqrt_ab'}
                    beta_t=sqrt(abs(alpha_s).*abs(beta)).*exp(i*0.5*(ph_beta+ph_alpha));

                case {'sqrt_abg'}
                    beta_t=(abs(alpha_s).*abs(beta).*abs(gamma_s)).^(1/3).*exp(i*(ph_beta+ph_alpha+ph_gamma)/3.0);

                case {'sqrt_ab_pa'}
                    beta_t=sqrt(abs(alpha_s).*abs(beta)).*exp(i*(ph_alpha));

                case {'sqrt_abg_pa'}
                    beta_t=(abs(alpha_s).*abs(beta).*abs(gamma_s)).^(1/3).*exp(i*ph_alpha);

                case {'max_ab'}
                    beta_t=max(abs(alpha_s),abs(beta)).*exp(i*(0.5*ph_beta+0.5*ph_alpha));

                case {'max_abg'}
                    beta_t=max(max(abs(alpha_s),abs(beta)),abs(gamma_s)).*exp(i*(ph_beta+ph_alpha+ph_gamma)/3.0);

                case {'max_ab_pa'}
                    beta_t=max(abs(alpha_s),abs(beta)).*exp(i*ph_alpha);

                case {'max_abg_pa'}
                    beta_t=max(max(abs(alpha_s),abs(beta)),abs(gamma_s)).*exp(i*ph_alpha);

                case {'min_ab_pa'}
                    beta_t=min(abs(alpha_s),abs(beta)).*exp(i*ph_alpha);    

                case {'avg_ab'}
                    beta_t=0.5*(alpha_s+beta);

                case {'avg_abg'}
                    beta_t=(1/3)*(alpha_s+beta+gamma_s);

                case {'avg_ab_pa'}
                    beta_t=0.5*(abs(alpha_s)+abs(beta)).*exp(i*(ph_alpha));

                case {'avg_abg_pa'}
                    beta_t=(1/3)*(abs(alpha_s)+abs(beta)+abs(gamma_s)).*exp(i*ph_alpha);

                case {'min_s_abg'}    
                    lb=[-1,-1,-1];
                    ub=[1,1,1];
                    options = optimset('Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-6);
                    norm=sum(abs(beta(:)).^4);

                    f=@(x)objfun_sharpness(x,abs(alpha_s),abs(beta),abs(gamma_s),norm);
                    [x fval]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
                    x=x/sum(x(:));
                    disp([x,fval])

                    beta_t=x(1)*abs(alpha_s)+x(2)*abs(beta)+x(3)*abs(gamma_s);
                    beta_t=beta_t.*exp(i*ph_alpha);

                case {'min_c_abg'}    
                    lb=[-1,-1,-1];
                    ub=[1,1,1];
                    options = optimset('Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-6);

                    f=@(x)objfun_chi(x,(alpha_s),(beta),(gamma_s),data);
                    [x fval]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
                    x=x/sum(x(:));
                    disp([x,fval])

                    beta_t=x(1)*alpha_s+x(2)*beta+x(3)*gamma_s;

                case {'min_TV_abg'}    
                    lb=[-1,-1,-1];
                    ub=[1,1,1];
                    options = optimset('Display','off','Algorithm','interior-point','TolFun',1e-6,'TolCon',1e-6);

                    f=@(x)objfun_TV(x,(alpha_s),(beta),(gamma_s));
                    [x fval]=fmincon(f,x0,[],[],[],[],lb,ub,[],options);
                    x=x/sum(x(:));
                    disp([x,fval])

                    beta_t=x(1)*alpha_s+x(2)*beta+x(3)*gamma_s;

            end


            if dims == 3, pn_g(:,:,qq)=beta_t;end
            if dims == 4, pn_g(:,:,:,qq)=beta_t;end

            beta=[];
            beta_t=[];
            ph_beta=[];
            alpha_s=[];

        end

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

c1=cross_correlation(a,conj_reflect(b));
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

