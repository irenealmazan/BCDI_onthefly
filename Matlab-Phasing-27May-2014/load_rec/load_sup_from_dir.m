function [ sup ] = load_sup_from_dir(ddir,flip)
%jclark
%loads a reconstruction
%from a parent directory only

if exist('flip') == 0,flip=0;end

%get sup fil
ampfile=rdir([ddir,'*SUP.rec']);

disp('Loading file ....')
disp(ampfile.name(numel(ddir)+1:end-4))

load(char(ampfile.name),'-mat')
sup=array;
array=0;

if flip ==1,sup=flipdim(sup,1);end


end



