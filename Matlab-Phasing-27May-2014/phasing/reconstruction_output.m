function reconstruction_output(filename, input_args )
%UNTITLED4 Summary of this function goes here
%  saves parameters from phasing to text file

fileid=fopen(filename,'w');

fields=(fieldnames(input_args));
nn=numel(fields);

input_args.coh='NA';

for ee=1:nn
    
    fprintf(fileid,'%s','--------------------');
    fprintf(fileid,'\t\n%s\t\n',char(fields(ee)));
    temp=getfield(input_args,char(fields(ee)));
    
    if iscellstr(temp) == 1
        for rr=1:numel(temp)
            fprintf(fileid,'\t%s\n',char(temp(rr)));
        end
    else
        if ischar(temp) == 1
            fprintf(fileid,'\t%s\n',(temp(:)));
        else
            fprintf(fileid,'\t%f\n',(temp(:)));
        end
    end
    
    
    
end    
    
fclose(fileid);

end

