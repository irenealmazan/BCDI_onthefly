function create_batch_phasingV2
%jclark
name_of_this_file='create_batch_phasingV2';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
%%

%only search for these strings to create the phaing for
big_strings={'*341727','*341849','*341924','*341947','*341957','*342038','*342066','*342137','*342148*','*342162*'};%,'*342169*'};

%how many files to break the phasing into
nfiles=3;

%deteremine how mnay per file
thpf=ceil(numel(big_strings)/nfiles);
lspf=nfiles*thpf-numel(big_strings);


%%
ind1=(1:thpf:numel(big_strings));
ind2=ind1+thpf-1;
ind2(end)=numel(big_strings);

for tt=1:nfiles
    
    strings=big_strings(ind1(tt):ind2(tt));

    fid = fopen([dir,'batch_phasing',num2str(tt),'.m'], 'w');

    for pp=1:numel(strings)

        script=['Matlabphasing_ver1_1.m'];
        n_sc=numel(script);
        files=rdir([[dir],'**/',char(strings(pp)),'/',script]);
        fnames=char(files.name);

        nn=size(fnames);



        for qq = 1:nn(1)

            sc_name=strtrim(fnames(qq,:));
            n_nm=numel(sc_name);


            [cur_dir file_out] = extract_dir_from_string(sc_name);

            fprintf(fid,'try \n');
            fprintf(fid,'\n');

            fprintf(fid,'disp(''Performing batch phasing...'') \n');

            fprintf(fid,'\n');
            fprintf(fid, 'cd ');
            fprintf(fid, cur_dir);
            fprintf(fid, '\n');

            fprintf(fid,'run ');
            fprintf(fid, cur_dir);
            fprintf(fid,script);
            fprintf(fid,'\n');
            cd(cur_dir);

            fprintf(fid,'disp(''complete...'') \n');
            fprintf(fid,'close all \n');
            fprintf(fid,'system(''purge''); \n');
            fprintf(fid,'end \n');
            fprintf(fid,'\n');


        end

    end


    fclose(fid);


end



end

