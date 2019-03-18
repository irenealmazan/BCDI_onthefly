function [ data ] =load_rec_data(params,save_dir,top_dir)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

try
    save_dir;
catch
    save_dir=[];
end
try
    top_dir;
catch
    top_dir=[];
end
try
    aliens=params.aliens;
catch
    aliens=0;
end
try
    nnc=params.nnc;
catch
    nnc=1;
end

files=params.files;
back=params.back;

bin=params.binning;
min_data=params.min_data;

file_params.no_hist=1;

if numel(top_dir) ~= 0

    for qq = 1:numel(files)
        
        name_t=char(files(qq));
        name_t=[name_t(1:end-3),lower(name_t(end-2:end))];
        name=rdir([top_dir,'**/*',name_t]);
        
        if numel(name) == 0,
            name_t=[name_t(1:end-3),upper(name_t(end-2:end))];
            name=rdir([top_dir,'**/*',name_t]);
        end
        
        full_files(qq)={name(1).name}; 
    end
           
    if numel(back) ~= 0,
        for qq = 1:numel(back)
            
            name_t=char(back(qq));
            name_t=[name_t(1:end-3),lower(name_t(end-2:end))];
            name=rdir([top_dir,'**/*',name_t]);
            
            if numel(name) == 0,
                name_t=[name_t(1:end-3),upper(name_t(end-2:end))];
                name=rdir([top_dir,'**/*',name_t]);
            end
            
            full_bg(qq)={name(1).name}; 
        
        end
    else full_bg=[];end

else
    data_dir=params.data_dir;

    full_files=strcat(data_dir,files);
    if numel(back) ~= 0,full_bg=strcat(data_dir,back);else full_bg=[];end
    
    end

data=bin_crop_center(full_files,full_bg,bin,min_data,aliens,nnc,params); %file_params

if numel(save_dir) ~= 0, save(save_dir,'data');end


end

