function savevtk1(array,params)

    params.ascii=1;

    fwid = fopen(params.vtkname,'w','b'); % IMPORTANT: big endian
    count = fprintf(fwid,'# vtk DataFile Version 2.0\n');
    count = fprintf(fwid,[params.comments,'\n']);
    if params.ascii
        count = fprintf(fwid,'ASCII\n');
    else
        count = fprintf(fwid,'BINARY\n');
    end
    count = fprintf(fwid,'DATASET STRUCTURED_POINTS\n');
    count = fprintf(fwid,'DIMENSIONS %u %u %u\n',params.dim);
    count = fprintf(fwid,'ORIGIN %u %u %u\n',params.origin);
    count = fprintf(fwid,'SPACING %3.2f %3.2f %3.2f\n',params.spacing);
    count = fprintf(fwid,'POINT_DATA %u\n',np);
    count = fprintf(fwid,['SCALARS ',params.varname,' ',vtkprecision,' 1\n']);
    count = fprintf(fwid,'LOOKUP_TABLE default\n');

    % write data to vtk file
    tic
    fprintf(1,['writing ',params.vtkname,' ... ']);
    if ~set.ascii
        count = fwrite(fwid, M, params.precision);
    else
        fprintf(fwid, '%g \n', M);
    end
    fclose(fwid);
    fprintf(1,'done in %5.3f s\n',toc);
end

