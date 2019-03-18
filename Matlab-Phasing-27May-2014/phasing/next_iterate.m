function pn = next_iterate(pn,pnm,support,ALG,beta,params)
% - Jesse Clark, LCN, UCL October-November 2010
%   jesse.clark@ucl.ac.uk, jessenclark@gmail.com
%calculates the mext iterate given a previous iterate, current guess and
%support.  uses different algorithms depending on whats specified with ALG
%available ones are given below.  just add another in if you require it
%ER,HIO,DM,SF,ASR,RAAR,HPR

switch ALG
    case {'GPHIO','GPHIOlr'}

    %if strcmp(ALG,'GPHIO','GPHIOlr')

        phi_max=max(params.phase_range);
        phi_min=min(params.phase_range);

        ph=angle(pnm);      %phase after modulus constraint
        ph_m=params.model_phase;   %model phase

        supportph=( ph < (ph_m+phi_max) & ph > (ph_m+phi_min ) );

        supportgp=support.*double(supportph); %create new support from old one and ph

        pn=(pnm.*supportgp)+(1-supportgp).*(pn-beta*pnm);

    %end

    
    case {'ER','ERlr'}
        %if strcmp(ALG,'ER','ERlr') 
        pn=pnm.*support;
        %end

    
    case {'HIO','HIOlr'}
    %if strcmp(ALG,'HIO','HIOlr')
        pn=(pnm.*support)+(1-support).*(pn-beta*pnm);
    %end  
    
    case {'OSS','OSSlr'}
        
        if isfield(params,'oss_alphax') ~= 1,params.oss_alphax=linspace(params.sz(2),1/params.sz(2),params.iterations);end
        if isfield(params,'oss_alphay') ~= 1,params.oss_alphay=linspace(params.sz(1),1/params.sz(1),params.iterations);end
        
        alphax=params.oss_alphax(params.itno);
        
        alphay=params.oss_alphay(params.itno);
        
        sigx=params.sz(2)/2/pi/alphax;
        sigy=params.sz(1)/2/pi/alphay;
        
        gx=ceil(3*sigx);
        gy=ceil(3*sigy);
        if gx > round(params.sz(2)/2),gx=round((params.sz(2)/2));end
        if gy > round(params.sz(1)/2),gy=round((params.sz(1)/2));end
        
        if ndims(pnm) == 2,
            gg=gauss_2D(gx,gy,sigx,sigy);end
        
        if ndims(pnm) == 3,
            if isfield(params,'oss_alphaz') ~= 1,params.oss_alphaz=linspace(params.sz(3),1/params.sz(3),params.iterations);end
            alphaz=params.oss_alphaz(params.itno);
            sigz=params.sz(3)/2/pi/alphaz;
            gz=ceil(3*sigz);
            if gz > round(params.sz(3)/2),gz=round((params.sz(3)/2));end
            gg=gauss_3D(gx,gy,gz,sigx,sigy,sigz);end
        
        out_s=convnfft((pnm.*support)+(1-support).*(pn-beta*pnm),gg,'same');
        pn=(pnm.*support)+(1-support).*out_s;
      
    

    case {'DM','DMlr'}
    %if strcmp(ALG,'DM','DMlr')
        pn=pn+pnm-pn.*support;
    %end

    case {'DMr','DMrlr'}
    %if strcmp(ALG,'DMr','DMrlr')
        pn=pn+.5*(pnm-pn.*support);
    %    pn=pnm+.5*(pn-pn.*support);
    %end

    case {'SF','SFlr'}
    %if strcmp(ALG,'SF','SFlr')
       pn=(2*support-1).*pnm ;
    %end


    case {'ASR','ASRlr'}
    %if strcmp(ALG,'ASR','ASRlr')
        pn=0.5*( (2*support-1).*(2*pnm-pn)+pn );
    %end
    
    case {'AAR','AARlr'}
    %if strcmp(ALG,'ASR','ASRlr')
        pn=0.5*( (2*support-1).*(2*pnm-pn)+pn );
    %end

    case {'RAAR','RAARlr'}
    %if strcmp(ALG,'RAAR','RAARlr')
        pn=0.5*beta*( (2*support-1).*(2*pnm-pn)+pn)+(1-beta)*pnm;
    %end

    case {'HPR','HPRlr'}
    %if strcmp(ALG,'HPR','HPRlr')
        pn=0.5*( (2*support-1).*(2*pnm-pn+(beta-1)*pnm)+pn+(1-beta)*pnm  );
    %end

    case {'HIO-OR','HIO-ORlr'}
    %if strcmp(ALG,'HIO-OR','HIO-ORlr')
        lama=.5;
        pn=(1+beta*(lama-1))*pn+(beta-lama-beta*lama)*support.*pn-beta*lama*pnm+(1+beta)*lama*pnm.*support;    
    %end 

    case {'GRAAR','GRAARlr'}
    %if strcmp(ALG,'GRAAR','GRAARlr')
        pn=0.5*beta*( (2*support-1).*(2*pnm-pn+(pn-pnm))+pn)+(1-beta)*pnm;
    %end

    case {'GHIO','GHIOlr'}
    %if strcmp(ALG,'GHIO','GHIOlr')
        pn=(pnm.*support)+(1-support).*(pn-beta*pnm)+.3*support.*(pn-pnm);
    %end 

    case {'HIO-ROR','HIO-RORlr'}
    %if strcmp(ALG,'HIO-ROR','HIO-RORlr')
        lama=random('uniform',-.5,.5);
        pn=(1+beta*(lama-1))*pn+(beta-lama-beta*lama)*support.*pn-beta*lama*pnm+(1+beta)*lama*pnm.*support;

    %end 


    case {'MEM','MEMlr'}
    %if strcmp(ALG,'MEM','MEMlr')
       if params.itno ~= 1,

           dpn=abs((pnm-pn)).*support;
           epsilon=1/max(abs(dpn(:)));
           %pn=beta*pnm.*exp(-dpn*epsilon)+(1-beta)*pnm;
           %pn=support.*pnm.*(.5*exp(-dpn*epsilon)+(1-.5));
           %pn=pnm.*cos(-dpn*epsilon).*support;
           pn=pnm.*exp(-dpn*epsilon).*support; 
       else
           pn=pnm.*support; 
       end
    %end


    %other algs, not conventional and (mostly) not published
    case {'SF-h-ER','SF-h-ERlr'}
    %if strcmp(ALG,'SF-h-ER','SF-h-ERlr')

       if params.itno < params.iterations/2
           disp('SF')
           pn=(2*support-1).*pnm ;    
       else
           disp('ER')
           pn=support.*pnm ;  
       end

    %end

    case {'ER-h-SF','ER-h-SFlr'}
    %if strcmp(ALG,'ER-h-SF','ER-h-SFlr')

       if params.itno > params.iterations/2
           disp('SF')
           pn=(2*support-1).*pnm ;    
       else
           disp('ER')
           pn=support.*pnm ;  
       end

    %end
    case {'ER-SF-ER','ER-SF-ERlr'}
   % if strcmp(ALG,'ER-SF-ER','ER-SF-ERlr')

       if params.itno > .1*params.iterations & params.itno < (params.iterations-5)
           disp('SF')
           pn=(2*support-1).*pnm ;    
       else
           disp('ER')
           pn=support.*pnm ;  
       end

    %end

    case {'ERs','ERslr'}
    %if strcmp(ALG,'ERs','ERslr') 
        pn=amp_sqrt(pnm.*support,0.5);
    %end

    case {'ER-hist','ER-histlr'}
    %if strcmp(ALG,'ER-hist','ER-histlr') 
        pn=histogram_constraint_v3(pn,pnm,support,params);
    %end

    case {'pcj-HIO','pcj-HIOlr'}
    %if strcmp(ALG,'pcj-HIO','pcj-HIOlr') 
        pnmpc=support.*phase_projector(pnm,params);
        pnpc=support.*phase_projector(pn,params);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end

    case {'HIO-hist','HIO-histlr'}
    %if strcmp(ALG,'HIO-hist','HIO-histlr')
        pnmpc=histogram_constraint_v3(pnm,pnm,support,params);
        pnpc=histogram_constraint_v3(pn,pn,support,params);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end

    case {'HIO-hist-sq','HIO-hist-sqlr'}
    %if strcmp(ALG,'HIO-hist-sq','HIO-hist-sqlr')
        pnmpc=hist_squared(pnm,support);
        pnpc=hist_squared(pn,support);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end

    case {'HIO-AMP','HIO-AMPlr'}
    %if strcmp(ALG,'HIO-AMP','HIO-AMPlr')
        pnmpc=amplitude_constraint(pnm,abs(pn));
        pnpc=amplitude_constraint(pn,abs(pn));    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end

    case {'ER-AMP','ER-AMPlr'}
    %if strcmp(ALG,'ER-AMP','ER-AMPlr')
        pnmpc=amplitude_constraint(pnm,abs(pn));    
        pn=pnmpc;
    %end

    case {'ERSF','ERSFlr'}
    %if strcmp(ALG,'ERSF','ERSFlr')
        if params.itno == 1,
            pn=support.*pnm;
        else
            pn=(2*support-1).*pnm ;
        end
    %end

    case {'SFa','SFalr'}
    %if strcmp(ALG,'SFa','SFalr')
       amp=(sqrt(abs(2*support.*pnm))-sqrt(abs(pnm))).^2;

       pn=(2*support-1).*pnm ;
       pn=amp.*atan2(imag(pn),real(pn));
    %end

    case {'HIOb','HIOblr'}
    %if strcmp(ALG,'HIOb','HIOblr')
        pnmpc=amp_boost(pnm.*support);
        pnpc=amp_boost(pn.*support);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end  

    case{'HIOs','HIOslr'}
    %if strcmp(ALG,'HIOs','HIOslr')
        pnmpc=amp_sqrt(pnm.*support,0.5);
        pnpc=amp_sqrt(pn.*support,0.5);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end 

    %if strcmp(ALG,'HIO-phgrd','HIO-phgrdlr')
    case {'HIO-phgrd','HIO-phgrdlr'}
        pnmpc=phase_gradient_constraint(pnm.*support);
        pnpc=phase_gradient_constraint(pn.*support);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end

    case {'HIOp','HIOplr'}
    %if strcmp(ALG,'HIOp','HIOplr')
        pnmpc=amp_sqrt_perc(pnm.*support,0.5,.02);
        pnpc=amp_sqrt_perc(pn.*support,0.5,.02);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end


    case {'HIOv','HIOvlr'}
    %if strcmp(ALG,'HIOv','HIOvlr')
        pnmpc=amp_sqrt_val(pnm.*support,0.5,.75);
        pnpc=amp_sqrt_val(pn.*support,0.5,.75);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end
    
    case {'RAARv','RAARvlr'}
    %if strcmp(ALG,'RAARv','RAARvlr')
        %pn=0.5*beta*( (2*support-1).*(2*pnm-pn)+pn)+(1-beta)*pnm;    
        %pnsm=amp_sqrt_val(pnm.*support,0.5,.7);
        %pns=amp_sqrt_val(pn.*support,0.5,.99);    
        pnsm=grey_level_my_image(pnm.*support,6);
        pns=grey_level_my_image(pn.*support,6);   

        pn=0.5*beta*(4*pnsm-2*pns-2*pnm+2*pn)+(1-beta)*pnm;

    %end

    case {'HIOso','HIOsolr'}
    %if strcmp(ALG,'HIOso','HIOsolr')

        if mod(params.itno,5) == 1
            pnmpc=amp_sqrt(pnm.*support,0.5);
            pnpc=amp_sqrt(pn.*support,0.5); 
        else
            pnmpc=pnm.*support;
            pnpc=pn.*support; 
        end
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end 

    case {'HIOsr','HIOsrlr'}
    %if strcmp(ALG,'HIOsr','HIOsrlr')
        pnmpc=amp_sqrt(pnm.*support,0.75);
        pnpc=amp_sqrt(pn.*support,0.75);    
        pn=pn+1*(pnmpc-beta*pnm-pnpc+beta.*pnmpc);
    %end 

    case {'HIOd','HIOdlr'}
    %if strcmp(ALG,'HIOd','HIOdlr')
        pnmpc=amp_sqrt_th(pnm.*support,0.5);
        pnpc=amp_sqrt_th(pn.*support,0.5);    
        pn=pnmpc+pn-beta*pnm-pnpc+beta.*pnmpc;
    %end 

    case {'HIOsi','HIOsilr'}
    %if strcmp(ALG,'HIOsi','HIOsilr')
        width=0.5;
        fact=params.itno/(width*params.iterations);
        if fact > 1,fact = 1;end

        if fact < 1,fact=0;end 

        disp(num2str(fact));
        pnmpc=fact*pnm.*support+(1-fact)*amp_sqrt(pnm.*support,0.5);
        pnpc=fact*pn.*support+(1-fact)*amp_sqrt(pn.*support,0.5);    
        pn=pn+.5*(pnmpc-beta*pnm-pnpc+beta.*pnmpc);
    %end 

end

end

