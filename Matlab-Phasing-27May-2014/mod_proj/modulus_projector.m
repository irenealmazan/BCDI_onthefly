function [pnm,error,params] = modulus_projector(pn,data,params,support)
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
%Modulus constraint projector
%Detailed explanation goes here
display=1;

try
    params.use_cf_mag;
catch
    params.use_cf_mag=0;
end

try
    params.angle;
catch
    params.angle=0;
    %if display == 1,disp(['DEFAULT - params.angle =',num2str(params.angle) ]),end
end

try
    params.use_fftconv;
catch
    params.use_fftconv=0;
    %if display == 1,disp(['DEFAULT - params.angle =',num2str(params.use_fftconv) ]),end
end
try
    params.mask;
catch
    params.mask=[];
end

try
    params.use_k_1;
catch
    params.use_k_1=0;
end

try
    params.use_2k_1;
catch
    params.use_2k_1=0;
end

try
    params.save_fft_pn;
catch
    params.save_fft_pn=0;
end

if params.use_2k_1 == 1,params.use_k_1=0;end  %2k-1 takes precedence over k-1

try
    params.mod_const_pow;
catch
    params.mod_const_pow=1;
end

try
    params.regularized_amp;
catch
    params.regularized_amp='none';
end
try
    params.silent;
catch
    params.silent=1;
end

try
    params.pcdi;
catch
    params.pcdi=0;
end

try
    params.no_zero;
catch
    params.no_zero=0;
end

%%
pnk=fftshift(fftn(fftshift(pn)));        %guess at wavefield in detector
Mk=abs(pnk);                             %modulus of iterate

if params.save_fft_pn == 1,params.Mk=Mk;end

Mm=abs(data);                            %measured modulus 

if params.pcdi == 1             %do PCDI related stuff
    
    if params.update == 1       %if it is update time do some things, otherwise don't to save time
        
        if params.pcdi_ER == 1,
            pnkER=fftshift(fftn(fftshift(pn.*support)));        %guess at wavefield in detector
            Mk=abs(pnkER);
            pnkER=0;            %modulus of iterate
        end

       if params.use_k_1 == 1,
            if params.silent ~= 1,disp('Using k-1 iterate....');end
            
            pnkER=fftshift(fftn(fftshift(params.pn)));        %guess at wavefield in detector
            Mk=abs(pnkER);
            pnkER=0;            %modulus of iterate
       end
       if params.use_2k_1 == 1,
            if params.silent ~= 1,disp('Using 2k-(k-1) iterate....');end
            pnkER=fftshift(fftn(fftshift(params.pn)));        %guess at wavefield in detector
            Mk=2*Mk-abs(pnkER);
            Mk(Mk<0)=0;
            pnkER=0;            %modulus of iterate
       end
    end
    
    %%%%USES deconvolution or model minimisation at each iteration to update the coh guess
   if params.update == 1,
       
       if params.silent ~= 1,disp('Updating coherence function....');end
       [coh,params]=determine_coh(Mk,Mm ,params);
   else
       if params.silent ~= 1 disp('Using previous coherence funciton....');end
       coh=params.coh;
   end
   
   Mk=abs(pnk);
   
   %11th Nov, changed to an fft baes gauss specific convolution functio in
   %gauss_kern_minimizer_ver2.  it is also changed below.
   
   switch lower(params.pcdi_type)
       case {'gauss','gauss_sa'}
       
         %if params.angle == 1000,Mk=sqrt(gauss_conv_fft(Mk.^2,params.kernalxy,params.angle));
         %else  Mk=sqrt(convn(Mk.^2,params.coh,'same'));end
         if params.use_fftconv == 0,
             if params.silent ~= 1,disp('Using direct convolution for modulus constraint....');end
             Mk=sqrt(convn(Mk.^2,params.coh,'same'));
         else
             if params.silent ~= 1,disp('Using FFT convolution for modulus constraint....');end
             Mk=sqrt(convnfft(Mk.^2,params.coh,'same'));
         end
         
       otherwise   
         %Mk=sqrt(convn(Mk.^2,params.coh,'same'));
         if params.use_fftconv == 0,
             if params.silent ~= 1,disp('Using direct convolution for modulus constraint....');end
             Mk=sqrt(convn(Mk.^2,params.coh,'same'));
         else
             if params.silent ~= 1,disp('Using FFT convolution for modulus constraint....');end
             Mk=sqrt(convnfft(Mk.^2,params.coh,'same'));
         end
   end
         
   sz=size(Mm);

   params.coh=coh;
   
 
end

if strcmp(lower(params.regularized_amp),'none')
    if numel(params.mask) <= 1
        
        %ratio=(pnk-pnk)+1;
        %ind=(Mk ~= 0);
        %ratio(ind)=(Mm(ind))./(Mk(ind));
        ratio = calc_ratio(Mm,Mk);  %calc the ratio, taking care of / by 0
        
        if params.no_zero == 1,
            ind=(Mm == 0);
            ratio(ind)=1;
        end

        if params.mod_const_pow ~= 1,      
            ratio=ratio.^(params.mod_const_pow);
        end

        pnkk=ratio.*pnk;   % aply modulus constraint  

    else

        if params.silent ~= 1
            disp('Using mask for saturated/missing pixels');end 

        
        if params.use_cf_mag == 0

            ratio = calc_ratio(Mm,Mk);
            ind=(params.mask ~= 0);
            ratio(ind)=1;
            pnkk=ratio.*pnk;   % aply modulus constraint  
    
        else
            
            pnkk=(2*Mm.^3-Mk.^3).*pnk;
            
            %pnkk=(2*Mm.-Mk.).*exp(i*angle(pnk));
            
        end
        
    %     amppnkk=abs(pnkk);
    %     phase=angle(pnkk);
    %     
    %     amp1=amppnkk.*(params.mask/max(params.mask(:)));
    %     amp1( amp1 <sqrt(max(params.mask(:))))=sqrt(max(params.mask(:)));
    %     amppnkk=params.mask/max(params.mask(:)).*amp1+(1-params.mask/max(params.mask(:))).*amppnkk;
    %     pnkk=amppnkk.*exp(i*phase);
    %     amp1=[];
    %     phase=[];
    %     ratio=[];   
    end
else
    
    Sq=data;            %gauss ~ pois ~ =sqrt(N)=sqrt(I)
    
    Eq=zeros(size(data));   %new func
    
    Sq( Sq ~= 0)=1/data(Sq ~= 0);  %1/s
    ratio = calc_ratio(Mm,Mk);
    
    if params.silent ~= 1
            disp(['Using regularized amplitude constraint -[ ',params.regularized_amp,' ]']);end 
    
    switch params.regularized_amp
        case 'uniform'
            Eq=1-1/24*(Sq).^2;
        case 'poisson'
            Eq=1-.236*(Sq).^(1.199);
        case 'gauss'
            Eq=1-1/8*(Sq).^2;
        otherwise
            Eq=1.0;
    end
    pnkk=Eq.*ratio.*pnk;
    Eq=[];
    Sq=[];
    ratio=[];
end

pnm=fftshift(ifftn(fftshift(pnkk)));        %return to sample plane
%%%%%%%%%%%%%
error=calc_chi(Mk(Mk ~= 0),Mm(Mk ~= 0));

end

function ratio = calc_ratio(Mm,Mk)

ratio=(Mm-Mm)+1;
ind=(Mk ~= 0);
ratio(ind)=(Mm(ind))./(Mk(ind));


end






