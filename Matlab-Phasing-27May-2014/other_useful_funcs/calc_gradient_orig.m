function [Ax,Ay,Az,Px,Py,Pz,ss] = calc_gradient_orig(amp,ph,do_s)
%jclark
%calcs gradient without wrapping issues for phase and amp.
%do_s is an optional imput to apply a support to theoutput to remove the
%step that arisrs if there is a boundry.  leave empty to not do or =0.  the
%value, 1,2,etc is how far from the boundry to remove in pixels.

if exist('do_s') ~= 1,do_s=0;end
if isempty('do_s'),do_s=0;end

%do amplitude first
aa=circshift(amp,[1,0,0]);
bb=circshift(amp,[0,0,0]);
Ay=aa-bb;

aa=circshift(amp,[0,1,0]);
bb=circshift(amp,[0,0,0]);
Ax=aa-bb;

if ndims(amp) == 3
    aa=circshift(amp,[0,0,1]);
    bb=circshift(amp,[0,0,0]);
    Az=aa-bb;

else Az=[];end

%do phase
Px=0;Pz=0;Py=0;
if exist('ph')
    
    if isempty('ph') ~= 1
        ph=exp(i*(ph));
        aa=circshift(ph,[1,0,0]);
        bb=circshift(ph,[0,0,0]);

        Py=angle(aa.*conj(bb)); %angle is equivalent to arg

        aa=circshift(ph,[0,1,0]);
        bb=circshift(ph,[0,0,0]);

        Px=angle(aa.*conj(bb));

        if ndims(ph) == 3
            aa=circshift(ph,[0,0,1]);
            bb=circshift(ph,[0,0,0]);

            Pz=angle(aa.*conj(bb));

        else Pz=[];end

    end
end

Px=-Px; %account for Jesse error

%do support
if do_s ~= 0
    support=shrink_wrap(abs(amp),.075,1.3);
    diffs=circshift(support,[do_s,0,0])-2*support;
    diffsy=(diffs == -1); 
    
    diffs=circshift(support,[0,do_s,0])-2*support;
    diffsx=(diffs == -1); 
     
    diffs=circshift(support,[0,0,do_s])-2*support;
    diffsz=(diffs == -1); 
    
    ss=diffsy.*diffsx.*diffsz;
end

    
end

