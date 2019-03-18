function  make_movies_mayavi_vtk(top_dir,frame_number,num_frames_per_second,final_size)
%jclark
%use for making a movie for the same view angle but different conditions,
%e.g. time, temperature, evaporation
%output movies from Make_moveAmp.py and phase

if exist('frame_number') ~= 1,frame_number=1;end
if isempty(frame_number) == 1,frame_number=1;end
if exist('final_size') ~= 1,final_size=500;end
if isempty(final_size) == 1,final_size=500;end
if exist('num_frames_per_second') ~= 1,num_frames_per_second =5;end
if isempty(num_frames_per_second) == 1,num_frames_per_second =5;end

ftype='.jpg';

numsize=3;      %how many places do the numbers got to 001 is 3 01 is 2 etc

base_dir=[top_dir,'movies/'];

if isdir(base_dir) == 0,mkdir(base_dir);end

kind={'Amp','Ph'};

for qq=1:numel(frame_number)
    for ww=1:numel(kind)
   
        numb=check_numlength(num2str(frame_number(qq)),numsize);
        
        string=[top_dir,'**/*',char(kind(ww)),'-',numb,ftype];
        save_dir=[base_dir,char(kind(ww)),'-',numb,'.avi'];
        fnames=sort_fnames_from_number(rdir(string),char(kind(ww)));
        make_avi_from_images(fnames,save_dir,num_frames_per_second,final_size);    
        
    end    
end



end

