function [data params] = load_MP_data(params)
%jclark
%returns the sqrt of the data
%will save the intensity though

try
    params.data_only;
catch
    params.data_only='NO';
end

try
    params.exit_on_error;
catch
    params.exit_on_error=1;
end

try
    params.binning;
catch
    params.binning=[1,1];
end

full_files=strcat(params.data_dir,params.files);
if numel(params.back) ~= 0,full_bg=strcat(params.data_dir,params.back);else full_bg=[];end

try
    data=bin_crop_center(full_files,full_bg,params.binning,params.min_data,params.aliens,params.nnc,params);
 catch
    
    if params.exit_on_error ~= 1
        disp(' ')
        disp('<<<<<<<<<<<<<<<<<<<<<<<ERROR>>>>>>>>>>>>>>>>>>>>>>>')
        disp('                  Error loading data               ')
        disp('<<<<<<<<<<<<<<<<<<<<<<<ERROR>>>>>>>>>>>>>>>>>>>>>>>')
        disp(' ')
        data=1;
    else
        error('**Error loading data.**  Check path and filenames are correct.  Set params.exit_on_error = 0 to continue without loading data.')
        
    end
end
    
params.nn=[size(data,2),size(data,1),size(data,3)];    %xyz

data=sqrt(abs(double(data)));

switch params.data_only
    
    case {'YES','YES-EXIT'}
        
        if isfield(params,'save_data_dir') ~= 1,params.save_data_dir=params.data_dir;end
        
        sname=[params.seq_let,'-',num2str(params.binning(1)),'X',num2str(params.binning(2))];
        disp(' ')
        disp('===========================================')
        disp('Saving prepared data (as the intensity)....')
        
        mat2tif(data.^2,[params.save_data_dir,'/Data-',sname,'.tif'])
        disp('Done....')
        disp('Saving as .vtk ....')
        savevtk2scalar(data.^2,[params.save_data_dir,'/Data-',sname,'.vtk'])
        disp('Done....')
        
        disp('Saving as .mat ....')
        array=data.^2;
        save([params.save_data_dir,'/Data-',sname,'.mat'],'array')
        array=[];
        disp('Done....')
        disp('===========================================')
end

if strcmp(params.data_only,'YES-EXIT')
    error('Exiting after data preparation.  Set params.data_only=''YES'' or ''NO'' to continue....')
end

end

