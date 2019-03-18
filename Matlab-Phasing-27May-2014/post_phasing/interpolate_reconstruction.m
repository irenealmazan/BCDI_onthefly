function [ output_args ] = interpolate_reconstruction(dir)
%UNTITLED Summary of this function goes here
%%   Detailed explanation goes here
try
    dir;
catch
    name_of_this_file='interpolate_reconstruction';
    dir_file=which(name_of_this_file);
    dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);
end
factor=2;
shrink=.5;
%%

fileAMP=rdir([[dir],'**/*AMP.rec']);
filePH=rdir([[dir],'**/*PH.rec']);
fileSUP=rdir([[dir],'**/*SUP.rec']);

AMP=load(fileAMP.name,'-mat');
AMP=AMP.array;
maxa=max(AMP(:));
mina=min(AMP(:));
AMP=zero_pad(AMP,size(AMP,2)*shrink,size(AMP,1)*shrink,size(AMP,3)*shrink);
AMP=FFT_interpolate(AMP,[factor*size(AMP,2),factor*size(AMP,1),factor*size(AMP,3)]);
AMP=(AMP-min(AMP(:)));
AMP=AMP/max(AMP(:))*(maxa-mina)+mina;

PH=load(filePH.name,'-mat');
PH=PH.array;
maxa=max(PH(:));
mina=min(PH(:));
PH=zero_pad(PH,size(PH,2)*shrink,size(PH,1)*shrink,size(PH,3)*shrink);
PH=FFT_interpolate(PH,[factor*size(PH,2),factor*size(PH,1),factor*size(PH,3)]);
PH=(PH-min(PH(:)));
PH=PH/max(PH(:))*(maxa-mina)+mina;

SUP=load(fileSUP.name,'-mat');
SUP=SUP.array;
SUP=zero_pad(SUP,size(SUP,2)*shrink,size(SUP,1)*shrink,size(SUP,3)*shrink);
SUP=FFT_interpolate(SUP,[factor*size(SUP,2),factor*size(SUP,1),factor*size(SUP,3)]);
SUP=SUP/max(SUP(:));
saveAMP=[fileAMP.name(numel(dir)+1:numel(fileAMP.name)-4),'-interp.rec'];
savePH=[filePH.name(numel(dir)+1:numel(filePH.name)-4),'-interp.rec'];
saveSUP=[fileSUP.name(numel(dir)+1:numel(fileSUP.name)-4),'-interp.rec'];

array=AMP;
save([dir,saveAMP],'array'); 

array=PH;
save([dir,savePH],'array'); 

array=SUP;
save([dir,saveSUP],'array'); 
array=0;

output_python_script(dir,'Mayavi_scene_script_Ver1-2-interp.py')

end

