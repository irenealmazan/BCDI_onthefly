function delete_random_phasing_script( top_dir )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

disp('Deleting Matlabphasing_ver1_1_Rnd#.m files....')
files=rdir([top_dir,'**/*','Matlabphasing_ver1_1_Rnd*.m']);

nfiles=max(size(files));

for qq=1:nfiles
    
    file=char(files(qq).name);
    disp('Deleting - ')
    disp(file)
    delete(file)
    disp(' ')
    
end

end

