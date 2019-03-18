function  spawn_random_reconstructions()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nstarts=5;

%get the directory
name_of_this_file='spawn_random_reconstructions';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

file='Matlabphasing_ver1_1_NOV_Batt2_AU.m';

rand_dir=[top_dir,'rand-starts/'];

%create a dir for the rand starts
if isdir(rand_dir) == 0,mkdir(rand_dir);end

%load the original script 
fid=fopen([top_dir,file]);
str=fscanf(fid,'%c');
fclose(fid);

%generate_names
seq_lets={'Rnd1'};
save_names={'Matlabphasing_ver1_1_Rnd1.m'};

for qq=2:nstarts,
    seq_lets=[seq_lets,['Rnd',num2str(qq)]];
    save_names=[save_names,['Matlabphasing_ver1_1_Rnd',num2str(qq),'.m']];
end

%check if there is already some files, if so delete these
prev=rdir([rand_dir,'*.m']);
if size(prev,1) > 0
   for qq=1:size(prev,1)
      delete(char(prev(qq).name))  
   end
    
end

%do the file stuff
for qq=1:nstarts
    
    old_str='name_of_this_file=''Matlabphasing_ver1_1''';
    new_str=['name_of_this_file=''Matlabphasing_ver1_1_Rnd',num2str(qq),''' '];
    new_file=regexprep(str,old_str,new_str,'once');
    
    old_str='data_dir=dir';
    new_str=['data_dir=''',top_dir,''''];
    new_file=regexprep(new_file,old_str,new_str,'once');
    
    old_str='seq_let=';
    new_str=['seq_let=','''Rnd',num2str(qq),''';%'];
    new_file=regexprep(new_file,old_str,new_str,'once');
    
    old_str='save_data=';
    new_str=['save_data=''NO'';%'];
    new_file=regexprep(new_file,old_str,new_str,'once');
    
    old_str='start_guess=';
    new_str='start_guess=''random'';%';
    new_file=regexprep(new_file,old_str,new_str,'once');
    
    
    fid=fopen([rand_dir,char(save_names(qq))],'w');
    fprintf(fid,'%c',new_file);
    fclose(fid);
    
end

fid = fopen([rand_dir,'random_batch_phasing.m'], 'w');
for qq = 1:nstarts
    
    fprintf(fid,'disp(''Performing batch phasing...'') \n');
    
    fprintf(fid,'\n');
    fprintf(fid, 'cd ');
    fprintf(fid, rand_dir);
    fprintf(fid, '\n');
    
    fprintf(fid,['rsd=',num2str(qq)]);
    fprintf(fid,'\ns = RandStream(''mcg16807'',''Seed'',rsd)');
    fprintf(fid,'\nRandStream.setDefaultStream(s)');
    
    fprintf(fid,'\nrun ');
    fprintf(fid, rand_dir);
    fprintf(fid,char(save_names(qq)));
    fprintf(fid,'\n');
    %cd(cur_dir);
    
    fprintf(fid,'disp(''complete...'') \n');
    %disp(sc_name)
    
end

cd(rand_dir)
run([rand_dir,'random_batch_phasing.m'])

end

