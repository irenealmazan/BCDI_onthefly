function [params] = calc_radial_stuff(amp,cryst)
%jclark
%amp is used for calculating the extent 
%cryst is the thing to calc as a funct of fractional size

nn=(size(amp));
mn=max(nn(:));

rad=1:mn;

SS=shrink_wrap(abs(amp),.075,1.3);

nshrink=100;

shrink=(4:nshrink)/nshrink;
shells=[];
for qq=1:numel(shrink)


   stemp=shrink_support(SS,shrink(qq));
   [Fx,Fy,Fz]=gradient(stemp);

   shell=abs(Fx)+abs(Fy)+abs(Fz);
   shell=(shell > 0);

   shells(:,:,:,qq)=shell;

   ind=find(shell == 1);

   nind(qq)=numel(ind);

   params.ud_mn_abs(qq)=mean(abs((cryst(ind))));

   params.ud_std_abs(qq)=std(abs((cryst(ind))));

   params.ud_mn(qq)=mean(((cryst(ind))));

   params.ud_std(qq)=std(((cryst(ind))));

   params.ud_rms(qq)=sqrt(mean(abs((cryst(ind))).^2));

end
        
params.shrink=shrink;

end

