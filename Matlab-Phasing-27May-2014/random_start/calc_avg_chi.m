function [avg] = calc_avg_chi(top_dir,search_string)
%used for calculating the average from random starts
% search_string is the search string
% eg. top_dir='/blah/blah/rand-starts/'
% search_string='**/*CVl*ERROR.mat'
% will return all error files from the reconstructions using CVl.
% search_string='**/*NM*ERROR.mat'
% for normal

files=rdir([top_dir,search_string]);


end

