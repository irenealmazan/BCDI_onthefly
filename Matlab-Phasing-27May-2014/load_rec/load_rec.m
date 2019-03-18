function pn = load_rec(amp_f,ph_f,flip)
%jclark
%loads a reconstruction given the amp and phase file names

try
    flip;
catch
    flip=0;
end


load(amp_f,'-mat')
amp=array;
array=0;
load(ph_f,'-mat');
ph=array;
array=0;

pn=amp.*exp(i*ph);

if flip == 1,pn=flipdim(pn,1);end


end

