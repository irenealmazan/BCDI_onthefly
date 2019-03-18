function copy_pngs_for_movie

string1='ROI_'; %

string2='X-Y'; %

ftype='.png';

auto_numb=1;        %give each one a number from 1 to N rather than getting it from a string

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_of_this_file='copy_pngs_for_movie';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
fdir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

fstring=[string1,string2,ftype];

fullnames=rdir([fdir,'**/*',fstring]);

nnames=size(fullnames,1);

save_dir=[fdir,'movies/',string1,string2,'/'];

if isdir(save_dir) ~= 1,mkdir(save_dir);end

for qq=1:nnames
    
    ffdir=extract_dir_from_string(char(fullnames(qq).name));
    
    if auto_numb ~= 1,numb=extract_number_from_string(ffdir);else numb=qq;end
    
    string=['cp ',char(fullnames(qq).name),' ',save_dir,num2str(numb),fstring];
    system(string)
    
end



end