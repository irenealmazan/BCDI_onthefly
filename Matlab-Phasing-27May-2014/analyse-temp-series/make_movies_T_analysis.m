function  make_movies_T_analysis(top_dir,num_frames_per_second,final_size)
%jclark
%output movies from tanalysis
if exist('final_size') ~= 1,final_size=500;end
if isempty(final_size) == 1,final_size=500;end
if exist('num_frames_per_second') ~= 1,num_frames_per_second =5;end
if isempty(num_frames_per_second) == 1,num_frames_per_second =5;end

ftype='.png';

base_dir=[top_dir,'movies/'];

if isdir(base_dir) == 0,mkdir(base_dir);end

orientation={'xy','xz','zy'};
kind={'I-Amp','I-Ph','S-Amp'};
save_prefix={'AMP','PH','SUP'};

for qq=1:numel(orientation)
    for ww=1:numel(kind)
   
        string=[top_dir,'**/*',char(kind(ww)),'-',char(orientation(qq)),ftype];
        save_dir=[base_dir,char(save_prefix(ww)),'-',char(orientation(qq)),'.avi'];
        make_avi_from_images(string,save_dir,num_frames_per_second,final_size);    
        
    end    
end



end

