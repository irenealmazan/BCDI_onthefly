name_of_this_file='conj_reflect_rec';

dir_file=which(name_of_this_file);
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);

fileAMP=rdir([[dir],'**/*AMP.rec']);

filePH=rdir([[dir],'**/*PH.rec']);

fileSUP=rdir([[dir],'**/*SUP.rec']);


AMP=load(fileAMP.name,'-mat');
AMP=AMP.array;
AMP=conj(ifftshift(fftn(fftshift(AMP))));
AMP=(ifftshift(ifftn(fftshift(AMP))));

PH=load(filePH.name,'-mat');
PH=PH.array;
PH=conj(ifftshift(fftn(fftshift(PH))));
PH=(ifftshift(ifftn(fftshift(PH))));


SUP=load(fileSUP.name,'-mat');
SUP=SUP.array;
SUP=conj(ifftshift(fftn(fftshift(SUP))));
SUP=(ifftshift(ifftn(fftshift(SUP))));


saveAMP=[fileAMP.name(numel(dir)+1:numel(fileAMP.name)-4),'-cjrf.rec'];
savePH=[filePH.name(numel(dir)+1:numel(filePH.name)-4),'-cjrf.rec'];
saveSUP=[fileSUP.name(numel(dir)+1:numel(fileSUP.name)-4),'-cjrf.rec'];

array=AMP;
save([dir,saveAMP],'array'); 

array=PH;
save([dir,savePH],'array'); 

array=SUP;
save([dir,saveSUP],'array'); 



