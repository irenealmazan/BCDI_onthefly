function savevtk2scalar(array,filename,array2,spacing)
%jclark
%adapted from savevtk.m, no does two scalars and set spacing
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    if exist('spacing') ~= 1,spacing=1;end
    
    disp('Writing VTK file....')
    disp(['Spacing - ',num2str(spacing)])
    
    %[nx, ny, nz] = size(array);
    [ny, nx, nz] = size(array);
   
    
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    fprintf(fid, 'DATASET STRUCTURED_POINTS\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid, '\n');
    fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
    fprintf(fid, 'SPACING    %d   %d   %d\n',spacing,spacing,spacing);
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS amp double\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');


    
    
    for a=1:nz
        
            fprintf(fid, '%d ', permute(array(:,:,a),[2 1]));
            %fwrite(fid,permute(array(:,:,a),[2 1]),'double');
            fprintf(fid, '\n');
        
    end
    
    if exist('array2') == 1
        if isempty(array2) ~= 1
            
            disp('Writing second scalar array to VTK file....')
            
            fprintf(fid, 'FIELD FieldData 1 \n');
            fprintf(fid, 'phases 1 %d', nx*ny*nz);
            fprintf(fid, 'double \n');

            for a=1:nz

                    fprintf(fid, '%d ', permute(array2(:,:,a),[2 1]));
                    %fwrite(fid, permute(array2(:,:,a),[2 1]),'double');
                    fprintf(fid, '\n');

            end    
        end
    end

    
    fclose(fid);
return

