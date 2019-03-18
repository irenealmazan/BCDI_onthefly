function [app_str] = make_string_from_recursions(before_string,strings,after_string)
%jclark
%make a list of strings that has a constant before and after but changing
%string (cell array, allows different size ones)

n_strings=max(size(strings));

for qq=1:n_strings
    
   app_str{qq}=[before_string,strings{qq},after_string];
       

end


end

