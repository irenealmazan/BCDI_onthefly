function [ pn ] = load_rec_from_dir(ddir,flip)
%jclark
%loads a reconstruction
%from a parent directory only

if exist('flip') == 0,flip=0;end

%get amp fil
ampfile=rdir([ddir,'*AMP.rec']);


%get phase fil
phfile=rdir([ddir,'*PH.rec']);

disp('Loading file ....')
disp(ampfile.name(numel(ddir)+1:end-4))

pn = load_rec(char(ampfile.name),char(phfile.name),flip);


end

