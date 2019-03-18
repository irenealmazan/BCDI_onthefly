function [fnames IIX] = sort_fnames_from_number(fnames,string_prefix)
%jclark
%will sort the names from rdir output accorindg ti numbers in the string
%put in a prefix to get the correct number

n_files=size(fnames,1);

for rr=1:n_files,
    numbs_init(rr)=extract_number_from_string(char(fnames(rr).name),[],string_prefix);
end

[BB IIX]=sort(numbs_init);

fnames=fnames(IIX);

end

