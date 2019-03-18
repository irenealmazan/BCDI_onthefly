function [ spacing ] = get_spacing_vtk(file)
%jclark

[ str_n ] = load_script_as_text(file,'POINT_DATA');

ind=regexp(str_n,'SPACING');

new_s=str_n(ind:end);       %get the bit with the spacing

cell_s=regexp(new_s,' ','split');  %split accordin to whitespace

nc=max(size(cell_s));

qq=0;
tt=0;
while qq == 0
    
    tt=tt+1;
    spacing=str2num(char(cell_s(tt)));
    qq=numel(spacing);
    
    
end

end

