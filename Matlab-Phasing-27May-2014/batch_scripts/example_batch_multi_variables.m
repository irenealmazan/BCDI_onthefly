
% example of doing a batch change of variable in Matlabphasing###.m
% list the old variables and there values as strings then list them with
% with the new values (see below).  remeber though the strings have to match
% exactly.

%current values within the file
%str_old={'params.GA=0;','iterations= 1000;','trigger=[5,950];','params.generations=1;','params.population=1000;'};
str_old={'params.sigma='};
%new values
%str_new={'params.GA=1;','iterations= 100;','trigger=[5,95];','params.generations=10;','params.population=100;'};
str_new={'params.sigma=.75;%'};


%% OR get it from the file location. comment out one
name_of_this_file='example_batch_multi_variables';
dir_file=which(name_of_this_file);    %finds the current folder for the phasing
top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory

%%
for qq=1:max(size(str_old))
    oldstr=char(str_old(1,qq));
    newstr=char(str_new(1,qq));
    batch_change_variable(top_dir,script,oldstr,newstr);
end