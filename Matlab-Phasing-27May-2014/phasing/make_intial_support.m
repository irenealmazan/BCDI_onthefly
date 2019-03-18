function [pn support] = make_intial_support(data,ss,params)
%UNTITLED3 Summary of this function goes here

if numel(size(data)) == 3
    
    support_1 = ones(ss);
    support = zeros(params.nn(1),params.nn(2),params.nn(3));
    support(params.nn(1)/2-size(support_1,1)/2+1:params.nn(1)/2+size(support_1,1)/2,params.nn(2)/2-size(support_1,1)/2+1:params.nn(2)/2+size(support_1,1)/2, ...
params.nn(3)/2-size(support_1,3)/2+1:params.nn(3)/2+size(support_1,3)/2)=support_1;
    %support=zero_pad_ver2(ones(ss),params.nn(1),params.nn(2),params.nn(3) );    %create support 
else
    support_1 = ones(ss);
    support = zeros(params.nn(1),params.nn(2),params.nn(3));
    support=zero_pad_ver2(ones([ss(1),ss(2)]),params.nn(1),params.nn(2) );
end

%if strcmp(params.start_guess,'flat'), pn=support;end%
%if strcmp(params.start_guess,'random'), pn=support.*random('uniform',0,1,[params.nn(2),params.nn(1),params.nn(3)]);end%

norm=sum(data(:).^2);

switch params.start_guess
   
    case 'flat'
        pn=support;
        disp('Using support as initial guess....')
        
    case 'random'
        pn=support.*random('uniform',0,1,[params.nn(2),params.nn(1),params.nn(3)]);
        disp('Using random support as initial guess....')
        
    case 'random-data'
        temp=data.*exp(i*2*pi*random('uniform',0,1,[params.nn(2),params.nn(1),params.nn(3)]));
        temp=ifftshift(ifftn(fftshift(temp)));
        pn=support.*temp;
        clear temp
        disp('Using data with random phase as intial guess....')
         
    case 'auto-sq'
        auto=fftshift(ifftn(fftshift(data.^2)));
        pn=abs(auto).^.5;
        auto=[];
        
    otherwise
    
        amp_f=params.start_guess;    
        ph_f=[amp_f(1:end-7),'PH.rec'];        
        pn=load_rec(amp_f,ph_f,1);
        load([ph_f(1:end-6),'SUP.rec'],'-mat')
        support=array;
        array=[];
        disp('Using previous reconstruction....')
        %need to check alignment
        F1=ifftshift(fftn(fftshift(pn)));
        F1=F1*sqrt(sum(abs(data(:).^2))/sum(abs(F1(:)).^2));
        disp('Aligning with data....')
        [h k l]=register_3d_reconstruction(abs(data).^2,abs(F1).^2);
        pn=ifftshift(ifftn(fftshift(sub_pixel_shift(F1,h,k,l))));
        
        support=flipdim(support,1);
        if is_conj_ref(pn,support) == 1,support=conj_reflect(support);end
        support=round(sub_pixel_shift(support,h,k,l));

        
        
        disp('Done....')
end




end

