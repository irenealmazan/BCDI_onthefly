function create_batch_spawn_rand_recons
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%%
name_of_this_file='create_batch_spawn_rand_recons';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
script='spawn_random_reconstructions.m';
n_sc=numel(script);
files=rdir([[dir],'**/*',script]);
fnames=char(files.name);

nn=size(fnames);

fid = fopen([dir,'batch_spawn_rand_recons.m'], 'w');

for qq = 1:nn(1)
    
    sc_name=strtrim(fnames(qq,:));
    n_nm=numel(sc_name);
    cur_dir=sc_name(1:end-n_sc);
    
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
    %disp(sc_name)
    
end

fclose(fid);

cd(dir)
run([dir,'batch_spawn_rand_recons.m'])


end

