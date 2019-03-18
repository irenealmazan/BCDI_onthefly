%
% Filename: $RCSfile: image_read_sub_help.m,v $
%
% $Revision: 1.1 $  $Date: 2008/06/10 17:05:14 $
% $Author: bunk $
% $Tag: $
%
% Description:
% parameter help for sub-routines of image_read like cbfread, edfread,
% speread and fliread
%
% Note:
% none
%
% Dependencies:
% none
%
%
% history:
%
% May 9th 2008, Oliver Bunk: 1st version
%
function [] = image_read_sub_help(m_file_name,extension,varargin)

% check minimum number of input arguments
if (nargin < 2)
    error('At least the m-file name and the extension have to be specified as input parameter.');
end

% accept cell array with name/value pairs as well
no_of_in_arg = nargin;
if (nargin == 3)
    if (isempty(varargin))
        % ignore empty cell array
        no_of_in_arg = no_of_in_arg -1;
    else
        if (iscell(varargin{1}))
            % use a filled one given as first and only variable parameter
            varargin = varargin{1};
            no_of_in_arg = 2 + length(varargin);
        end
    end
end

% check number of input arguments
if (rem(no_of_in_arg,2) ~= 0)
    error('The optional parameters have to be specified as ''name'',''value'' pairs');
end


% parse the variable input arguments
examples = 1;
extension_returned = 0;
for ind = 1:2:length(varargin)
    name = varargin{ind};
    value = varargin{ind+1};
    switch name
        case 'Examples' 
            examples = value;
        case 'ExtensionReturned'
            extension_returned = 1;
        otherwise
            error('unknown parameter %s',name);
    end
end

fprintf('Usage:\n')
if (extension_returned)
    fprintf('[frame]=%s(<filename> [[,<name>,<value>] ...]);\n',...
        m_file_name);
else
    fprintf('[frame]=%s(<filename> [[,<name>,<value>] ...]);\n',...
        m_file_name);
end
fprintf('The optional <name>,<value> pairs are:\n');
fprintf('''RetryReadSleep'',<seconds>           if greater than zero retry opening after this time (default: 0.0)\n');
fprintf('''RetryReadMax'',<0-...>               maximum no. of retries, 0 for infinity (default: 0)\n');
fprintf('''MessageIfNotFound'',<0-no,1-yes>     display a mesage if not found, 1-yes is default\n');
fprintf('''ErrorIfNotFound'',<0-no,1-yes>       exit with an error if not found, default is 1-yes\n');
if (examples)
    fprintf('\n');
    fprintf('Examples:\n');
    fprintf('[frame]=%s(''~/Data10/roper/image.%s'');\n',...
        m_file_name,extension);
    fprintf('\n');
    fprintf('The returned structure has the fields data, header and extension.\n');
end
