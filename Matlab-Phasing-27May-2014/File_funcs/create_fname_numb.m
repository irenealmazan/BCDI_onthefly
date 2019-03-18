function names=create_fname_numb(prefix,suffix,numbers)

nfiles=max(size(numbers));

names='';

nzeros=numel(num2str(numbers(end)));


for qq = 1:nfiles
    
    nnzro='';
    while (numel(nnzro)+numel(num2str(numbers(qq)))) < nzeros
        nnzro=nnzro+'0';
    end
    
    name=[prefix,nnzro,num2str(numbers(qq)),suffix];
    names=[names,{name}];
    
end




end

