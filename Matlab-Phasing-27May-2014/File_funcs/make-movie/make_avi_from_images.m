function make_avi_from_images(string,save_dir,num_frames_per_second,final_size,N)
%jclark
%writes an avi from a string specifiying a recusion relation or cell
%containing files
%pass through a string to search for those files or alternatively pass
%through file names using a cell array
%this avi writer works on the mac

%[filename, pathname] = uigetfile('*.*','MultiSelect', 'on','Pick a file');

if exist('final_size') ~= 1,final_size=500;end
if isempty(final_size) == 1,final_size=500;end
if exist('num_frames_per_second') ~= 1,num_frames_per_second =5;end
if isempty(num_frames_per_second) == 1,num_frames_per_second =5;end


%get the files
if isstruct(string) ~= 1
    fnames=rdir([string]);
else 
    fnames=string;
end

if exist('N') ~= 1,N=size(fnames,1);end
if isempty(N) == 1,N=size(fnames,1);end


if N ~= 0
    %Precharge some variables to improve speed
    name1=[char(fnames(1).name)];
    [a,map]=imread(name1);
    %mov(1:N) = struct('cdata',a,'colormap',map);

    sz=size(a);
    frac=final_size/max(sz(:));

    %aviobj = avifile ( [save_dir], 'fps',num_frames_per_second); 
    myVideo = VideoWriter(save_dir);
    myVideo.FrameRate = num_frames_per_second;
    myVideo.Quality = 100;
    open(myVideo);

    %Main loop to convert images to frames and add them to movie.
    for qq=1:N
        name1=[char(fnames(qq).name)];
        [a,map]=imread(name1);
        a=imresize(a,frac);

        %mov(qq)=im2frame(a,map);
        %aviobj = addframe ( aviobj,a);

        writeVideo(myVideo,a);
    end
    %Convert movie to avi file and save it on current folder
    %movie2avi(mov, [save_dir,'.avi'], 'compression', 'None','fps',5);
    %aviobj = close ( aviobj );
    close(myVideo);
    
else
    disp(' ')
    disp('Could not find files....')
    disp(' ')
end

end