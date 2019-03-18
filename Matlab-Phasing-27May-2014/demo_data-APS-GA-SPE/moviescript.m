%% movie of the timing scan
cd ~/Desktop/
    vidObj = VideoWriter('pdafterdefect.avi');
    open(vidObj);
 
% Create an animation.
pos = get(gcf,'position');
 
nm=1;
 for i=1:size(data,3)
        imagesc(log10(data(25:225,25:225,i,nm))); axis equal; title(['Angle number ' num2str(round(i))]); axis off;
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
        writeVideo(vidObj,currFrame);
         writeVideo(vidObj,currFrame);
          writeVideo(vidObj,currFrame);
 end
 
   close(vidObj);