function [ data ] = load_ND_data(params)
%jclark
%to load data with the option of 4d (3d ptycho)


try
    params.load_seperately; %check if data sets should be combined or loaded seperately (4d -> 3d ptycho)
catch
    params.load_seperately=0;
end

%load_seperately=0 will just sum data sets
if params.load_seperately == 0
    full_files=strcat(params.data_dir,params.files);
    if numel(params.back) ~= 0,full_bg=strcat(params.data_dir,params.back);else full_bg=[];end
    data=bin_crop_center(full_files,full_bg,params.bin,params.min_data,params.aliens,params.nnc,params);
else
    ndata=size(params.files,2);
    
    for qq = 1:ndata
    
        full_files=strcat(params.data_dir,params.files(qq));
        if numel(params.back) ~= 0,full_bg=strcat(params.data_dir,params.back(qq));else full_bg=[];end
        adata=bin_crop_center(full_files,full_bg,params.bin,params.min_data,params.aliens,params.nnc,params);
        
        if qq == 1,data=zeros([size(adata),ndata]);end
        
        switch ndims(data)
            case 4
                data(:,:,:,qq)=adata;
                adata=[];
            case 3
                data(:,:,qq)=adata;
                adata=[];
        end
    end
    
end



end

