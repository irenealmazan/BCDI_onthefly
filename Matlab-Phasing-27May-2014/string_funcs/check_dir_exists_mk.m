function check_dir_exists_mk( ddir)
%jclark

if isdir(extract_dir_from_string(ddir)) ~= 1,
    
    mkdir(extract_dir_from_string(ddir))
    
end

end

