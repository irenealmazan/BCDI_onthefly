function copy_multi_pngs_for_movie_rotating

strings1={'Ph-zy','Ph-zy','Ph-zy'};
strings2={'-24','-18','-12'};

nstrings=size(strings1,2);

for tt=1:nstrings
    
    string1=strings1{tt}; %

    string2=strings2{tt}; %

    ftype='.png';

    nlen=3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    name_of_this_file='copy_multi_pngs_for_movie_rotating';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    fdir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

    fstring=[string1,string2,ftype];

    fullnames=rdir([fdir,'**/*',fstring]);

    nnames=size(fullnames,1);

    save_dir=[fdir,'movies/',string1,string2,'/'];

    if isdir(save_dir) ~= 1,mkdir(save_dir);end

    for qq=1:nnames

        ffdir=extract_dir_from_string(char(fullnames(qq).name));
        numb=extract_number_from_string(ffdir);
        numb=check_numlength(numb,nlen);
        string=['cp ',char(fullnames(qq).name),' ',save_dir,(numb),fstring];
        system(string)

    end

end


end