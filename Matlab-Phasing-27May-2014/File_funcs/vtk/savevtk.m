function savevtk(array, filename)
%  savevtk Save a 3-D scalar array in VTK format.
%  savevtk(array, filename) saves a 3-D array of any size to
%  filename in VTK format.
    
    [ny, nx, nz] = size(array);
    dataType='float';
    x=(1:nx)-nx/2;
    y=(1:ny)-ny/2;
    z=(1:nz)-nz/2;
    
    fid = fopen(filename, 'wt');
    fprintf(fid, '# vtk DataFile Version 2.0\n');
    fprintf(fid, 'Comment goes here\n');
    fprintf(fid, 'ASCII\n');
    fprintf(fid, '\n');
    
    %fprintf(fid, 'DATASET STRUCTURED_ POINTS\n');
    %fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    %fprintf(fid, '\n');
    %fprintf(fid, 'ORIGIN    0.000   0.000   0.000\n');
    %fprintf(fid, 'SPACING    1.000   1.000   1.000\n');
    
    fprintf(fid, 'DATASET RECTILINEAR_ GRID\n');
    fprintf(fid, 'DIMENSIONS    %d   %d   %d\n', nx, ny, nz);
    fprintf(fid,'X_ COORDINATES %d   %s \n',nx,  dataType);
    fprintf(fid,'%g \n',x);
    fprintf(fid,'Y_ COORDINATES %d   %s \n',ny,  dataType);
    fprintf(fid,'%g \n',y);
    fprintf(fid,'Z_ COORDINATES %d   %s \n',nz,  dataType);
    fprintf(fid,'%g \n',z);
    
    fprintf(fid, '\n');
    fprintf(fid, 'POINT_DATA   %d\n', nx*ny*nz);
    fprintf(fid, 'SCALARS scalars double\n');
    fprintf(fid, 'LOOKUP_TABLE default\n');
    fprintf(fid, '\n');
    
    fprintf(fid, '%g \n', array);
    %for a=1:nz
    %    for b=1:ny
    %        for c=1:nx
    %            fprintf(fid, '%d ', array(c,b,a));
    %        end
    %        fprintf(fid, '\n');
    %    end
    %end
    
    fclose(fid);
return


