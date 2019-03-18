function [ kernalxy angle] = gauss_kern_minimizer_ver2(slice,modulus,kernalxy,angle,xnotdo,ynotdo,znotdo,anglenotdo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%def gauss_kern_minimizer(slice,modulus,kernalxy=[],xnotdo=[],ynotdo=[],roixy=[],xy_range=[],norm=[],things=[]):

%JESSE CLARK AUGUST 2010
%slice - the unconvolved estimate
%modulus -  the convolved data
%ynotdo,xnotdo - set to 1 will not do the respective axis
%roixy - the roi over which to compare the two, [+-x,+-y] default is [32,32]
%xy_Range - the range over which to search, given as [min,max]
%norm - renormalizes prior to comparing
%things - the number of steps, default is 11

%best usage is to do two or more instances of gauss_kern_minimizer with a shrinking roi
%eg.
%   kernalxy=gauss_kern_minimizer(data,modulus,xnotdo=xnotdo,ynotdo=ynotdo)%,things=6)%,kernalxy=kernalxy)
%   kernalxy=gauss_kern_minimizer(data,modulus,xnotdo=xnotdo,ynotdo=ynotdo,kernalxy=kernalxy,xy_range=[.2,.2])
%   kernalxy=gauss_kern_minimizer(data,modulus,xnotdo=xnotdo,ynotdo=ynotdo,kernalxy=kernalxy,xy_range=[.2,.2])
%this seems to work, although there is something not wuite right.  for some values if sig the summed diff
%is the same for two different sets of values.  i think its a preciosin problem but increasing it doesn't seem to work
    
i=complex(0.,1.);
sx=size(slice)  ;                                        %get size   

                
gn=11                    ;                               %size of the guass kernal for convolving

%anglenotdo=0    ;    
mod_comp=modulus;           %get roi, remember its row column major
absF2=slice;

details=1;

%totF2=sum(sum(sum(mod_comp.^2)));   

%% 3D version
if ndims(slice) == 3,
    angle=0;
    %these control the way the algorithm searches.  can playa round with these
    %if you want
    things=11;%7;%11             needs to be odd   ;                               %the number of points for the line search
    shrink=[1,5];%[1,4,20]; %[1,10.0]              ;                               %the factor by which to shrink in each progressive cycle
    search_range=.5;%0.6;  %.5                                     %sets the initial search range

    rounds=max(size(shrink))                    ;                               %the number of cycles to do the line search
    
    if details == 1, disp('Doing 3D minimisation....'), end
    params=3;

    for ww = 1:rounds

        min_range=kernalxy-search_range/shrink(ww);
        max_range=kernalxy+search_range/shrink(ww);

        sigs=zeros(params,things);

        if xnotdo ~=0, sigs(1,:)=0;end
        if ynotdo ~=0, sigs(2,:)=0;end
        if znotdo ~=0, sigs(3,:)=0;end

        for aa= 1:params,sigs(aa,:)=(0:(things-1))/(things-1).*(max_range(aa)-min_range(aa))+min_range(aa);end

        ind=((things+1.)/2.)   ;                         %the location of the min value for sig    

        % do x
        if xnotdo == 0
            if ww == 1, if details == 1, disp('Doing x minimisation....'), end,end
            sigx=min_gauss_param3d(gn,gn,gn,sigs(1,:),sigs(2,ind),sigs(3,ind),mod_comp,absF2);
        else sigx=0;end
        % do y
        if ynotdo == 0
            if ww == 1, if details == 1, disp('Doing y minimisation....'), end,end
            sigy=min_gauss_param3d(gn,gn,gn,sigx,sigs(2,:),sigs(3,ind),mod_comp,absF2);
        else sigy=0;end
        % do z
        if znotdo == 0
            if ww == 1, if details == 1, disp('Doing z minimisation....'), end,end
            sigz=min_gauss_param3d(gn,gn,gn,sigx,sigy,sigs(3,:),mod_comp,absF2);
        else sigz=0;end

        kernalxy=[sigx,sigy,sigz];

    end
end
%%
%% 2D version
if ndims(slice) == 2,
    if details == 1,disp('Doing 2D minimisation....'),end
    params=2;
    kernalxy=[kernalxy(1),kernalxy(2)];
    
    %these control the way the algorithm searches.  can playa round with these
    %if you want
    things=11;%             needs to be odd   ;                               %the number of points for the line search
    shrink=[1,5] ;             ;                               %the factor by which to shrink in each progressive cycle
    search_range=.5;                                     %sets the initial search range

    rounds=max(size(shrink))                    ;                               %the number of cycles to do the line search
    
    for ww = 1:rounds

        min_range=kernalxy-search_range/shrink(ww);
        max_range=kernalxy+search_range/shrink(ww);

        err=zeros(params,things)   ;                   %make an array for the error values

        sigs=zeros(params,things);
        
        max_ang=angle+20/shrink(ww);
        min_ang=angle-20/shrink(ww);
        angs=(0:(things-1))/(things-1).*(max_ang-min_ang)+min_ang;
        
        if xnotdo ~=0, sigs(1,:)=0;end
        if ynotdo ~=0, sigs(2,:)=0;end
        if anglenotdo ~=0, angs(:)=0;end

        for aa= 1:params,sigs(aa,:)=(0:(things-1))/(things-1).*(max_range(aa)-min_range(aa))+min_range(aa);end

        ind=((things+1.)/2.)   ;                         %the location of the min value for sig    

        % do x
        if xnotdo == 0
            if ww == 1, if details == 1, disp('Doing x minimisation....'), end,end
            sigx=min_gauss_param2d(gn,gn,sigs(1,:),sigs(2,ind),angle,mod_comp,absF2);
        else sigx=0;end
        % do y
        if ynotdo == 0
            if ww == 1, if details == 1, disp('Doing y minimisation....'), end,end
            sigy=min_gauss_param2d(gn,gn,sigx,sigs(2,:),angle,mod_comp,absF2);
        else sigy=0;end
        % do angle
        if anglenotdo == 0
            angle=min_gauss_param2d(gn,gn,sigx,sigy,angs,mod_comp,absF2);
        end
        
        sigz=0;
        kernalxy=[sigx,sigy,sigz];

    end
end


end



function [min_val]=min_gauss_param2d(nx,ny,sigx,sigy,angle,mod_comp,absF2)

sigx0=squeeze(sigx);
sigy0=squeeze(sigy);
angle0=squeeze(angle);

sx=max(max(size(sigx)));
sy=max(max(size(sigy)));
sa=max(max(size(angle0)));


loc=find([sx,sy,sa] == max([sx,sy,sa])); %find which one is the one to test values for

things=max([sx,sy,sa]);  %number to iterat over

if sx == 1, sigx=zeros([1,things])+sigx0;else sigx=sigx0;end    %if there is only one element create ana rray full of them
if sy == 1, sigy=zeros([1,things])+sigy0;else sigy=sigy0;end
if sa == 1, angle=zeros([1,things])+angle0;else angle=angle0;end

error=zeros(1,things);    

ind_z=find(mod_comp);                   %get nonzero values
totF2=sum(sum(sum(absF2.^2)));

for qq = 1:things
    
    kernal=gauss_2D(nx,ny,sigx(qq),sigy(qq),angle(qq));
    
    
    guess=convn(absF2.^2,kernal,'same');
    %guess=gauss_conv_fft(absF2.^2,[sigx(qq),sigy(qq)],angle(qq) );
    
    diff=abs(guess(ind_z)-mod_comp(ind_z).^2);     

    error(qq)=sum(sum(sum(diff)))/totF2   ;

end
   
ind=find(error == min(error));              %locate minimum 
ind=ind(1);
min_val=[sigx(ind),sigy(ind),angle(ind)];
min_val=min_val(loc);                       %return the correct one


end

function [min_val]=min_gauss_param3d(nx,ny,nz,sigx,sigy,sigz,mod_comp,absF2)

sigx0=squeeze(sigx);
sigy0=squeeze(sigy);
sigz0=squeeze(sigz);

sx=max(max(size(sigx)));
sy=max(max(size(sigy)));
sz=max(max(size(sigz)));

loc=find([sx,sy,sz] == max([sx,sy,sz])); %find which one is the one to test values for
things=max([sx,sy,sz]);  %number to iterat over

if sx == 1, sigx=zeros([1,things])+sigx0;else sigx=sigx0;end    %if there is only one element create ana rray full of them
if sy == 1, sigy=zeros([1,things])+sigy0;else sigy=sigy0;end
if sz == 1, sigz=zeros([1,things])+sigz0;else sigz=sigz0;end

error=zeros(1,things);    

ind_z=(mod_comp > 0 );                   %get nonzero values
totF2=sum(sum(sum(absF2.^2)));

for qq = 1:things
    
    kernal=gauss_3D(nx,ny,nz,sigx(qq),sigy(qq),sigz(qq));
    %guess=convn(absF2.^2,kernal,'same');
    guess=gauss_conv_fft(absF2.^2,[sigx(qq),sigy(qq),sigz(qq)]);
    diff=abs(guess(ind_z)-mod_comp(ind_z).^2);     

    error(qq)=sum(sum(sum(diff)))/totF2   ;

end
   
ind=find(error == min(error));              %locate minimum 
ind=ind(1);
min_val=[sigx(ind),sigy(ind),sigz(ind)];
min_val=min_val(loc);                       %return the correct one


end
