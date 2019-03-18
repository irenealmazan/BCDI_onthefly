function [ pn_g supports] = GA_pre_screen(support,data,params)
%jclark
npop=params.population;

if isfield(params,'GA_prescreen_n') ~= 1,
    params.GA_prescreen_n=10*npop;end  %how many to check

if numel(params.sz) == 3
   gg=gauss_3D(11,11,11,1,1,1); 
else
   gg=gauss_2D(11,11,1,1);
end

for qq=1:params.GA_prescreen_n

    if params.GA_data_rand == 0
        pni=support.*random('uniform',.95,1,[params.sz]);
        temp=fftshift(fftn(fftshift(pni)));
        temp=data.*exp(i*angle(temp));
        pni=fftshift(ifftn(fftshift(temp))); %make the first iterate
    else
        %temp=data.*exp(i*2*pi*random('uniform',0,1,[params.sz]));
        temp=data.*exp(i*2*pi*random('uniform',0,1,[params.sz]) );
        
        pni=support.*ifftshift(ifftn(fftshift(temp)));        
    end

    sharpness(qq)=sum(abs(pni(:)).^4); 
    
    if qq <= npop                     %fill the first with anything
        
        if numel(params.sz) == 3
            pn_g(:,:,:,qq)=pni;
            supports(:,:,:,qq)=support;
        else
            pn_g(:,:,qq)=pni;
            supports(:,:,qq)=support;
        end
        
        sharp_keep(qq)=sharpness(qq);
        
        
        
    else                              %only check if there is already some in there
        
        ind=find(sharp_keep == max(sharp_keep(:))); %find worst one
        
        if sharpness(qq) < sharp_keep(ind)       %is current one better than worst?
            
            if numel(params.sz) == 3
                pn_g(:,:,:,ind)=pni;        %replace worst one
                sharp_keep(ind)=sharpness(qq); 
            else
                pn_g(:,:,ind)=pni;
                sharp_keep(ind)=sharpness(qq); 
            end
            
        end
        
    end

    
end


if params.GA_prescreen_n > npop

end



end

