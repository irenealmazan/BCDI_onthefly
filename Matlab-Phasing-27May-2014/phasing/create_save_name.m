function [ name ] = create_save_name( files,ALG1,ALG2,iterations,pcdi,sw ,seq,type,GPU,params)
%jclark
%files are the files that were used
%ALG1 and ALG2 are the algorithms
%iterations -
%pcdi, weather pcdi correction is applied 1 yes 0 no
%sw, flag that tells weather sw was on at any stage
%seq identifier to do multiple reconstructions 
%type - type of pcdi correction

%get file number from the spe data file
try
    params.GA;
catch
    params.GA=0;
end

try
    params.do_2D;
catch
    params.do_2D=0;
end

try
    params.start_guess;
catch
    params.start_guess='flat';
end

try
    params.GA_metric;                 %how to decide what is the 'best'
catch
    params.GA_metric='chi';                 %how to decide what is the 'best'
end
try
    params.iterate_avg_enabled;
catch
    params.iterate_avg_enabled=0;
end
try
    params.regularized_amp;
catch
    params.regularized_amp='none';
end

%%
a=cell2mat(files(1));

try
    if isempty(strfind(lower(char(files{1})),'.tif')) == 0
        idx = regexp(char(files{1}),'\d+');
        nums = regexp(char(files{1}),'\d+','match');   %get the numbers from the file name
        numb=num2str(char(nums));
        szn=size(numb);
        if szn(1) ~= 1,numb=numb(end,:);end
        numb=strtrim(numb);
    end

    if isempty(strfind(lower(char(files{1})),'.spe')) == 0
        numb=num2str(sscanf(a(strfind(a,'-')+1:numel(a)),'%i'));
        szn=size(numb);
        if szn(1) > 1,numb=numb(end,:);end
        numb=strtrim(numb);
    end
    if isempty(strfind(lower(char(files{1})),'.mat')) == 0
        idx = regexp(char(files{1}),'\d+');
        nums = regexp(char(files{1}),'\d+','match');   %get the numbers from the file name
        numb=num2str(char(nums));
        szn=size(numb);
        if szn(1) ~= 1,numb=numb(end,:);end
        numb=strtrim(numb);
    end
catch
    numb='';
end

its=num2str(iterations);

try 
    GPU;
catch
    GPU=0;
end

if pcdi == 1, 
    
    switch lower(type)
        
        case 'lucy'
            pc='CVl';
        
        case 'lucy_ps'
            pc='CVlps';
            
        case 'lucy_lsqs'
            pc='CVlsq';
        
        case 'weiner'
            pc='CVw';
            
        case 'lucy-weiner'
            pc='CVlw';
            
        case 'gauss'
            pc='CVg';
        
        case 'gauss_ps'
            pc='CVgps';
            
        case 'gauss_sa'
            pc='CVgsa';
            
        case 'gauss_um'
            pc='CVgum';
            
        case 'gauss_cm'
            pc='CVgcm';
            
        case 'lsqs'
            pc='CVlsq';
            
        case 'HyBR'
            pc='CVHyBR';
            
        otherwise
            pc='CV?';
            
            
    end
    %if strcmp(type,'lucy'),pc='CVl';end
    %if strcmp(type,'gauss'),pc='CVg';end    
else pc='NM';end

if sw == 1, 
    son='SW';
    if strcmp(params.sw_type,'percent'), son='PC';end
    
    if strcmp(params.sw_type,'gauss_percent'), son='PC';end
    
else son='FS';end

% if GPU == 1,
%     name=['Rec-',seq,'-',numb,'-GPU-',ALG1,ALG2,its,'-',pc,'-',son];
% else name=['Rec-',seq,'-',numb,'-',ALG1,ALG2,its,'-',pc,'-',son];end

name=['Rec-',seq,'-',numb,'-'];

switch params.regularized_amp
    
    case 'gauss'
        ra='Rg';
        name=[name,ra,'-'];   
    case 'poisson'
        ra='Rp';
        name=[name,ra,'-'];   
    case 'uniform'
        ra='Ru';
        name=[name,ra,'-'];   
end

if GPU == 1
    name=[name,'GPU-']; 
end

if params.GA == 1
    
    switch params.breed_mode1
        
        case 'sqrt_ab'  
            GA_str='s';
        case 'sqrt_ab_pa'  
            GA_str='t';
        case 'sqrt_abg'  
            GA_str='u';
        case 'sqrt_abg_pa'
            GA_str='v';
        
        case 'max_ab'  
            GA_str='m';
        case 'max_ab_pa'  
            GA_str='n';
        case 'max_abg'  
            GA_str='o';
        case 'max_abg_pa'
            GA_str='p';
            
        case 'avg_ab'  
            GA_str='a';
        case 'avg_ab_pa'  
            GA_str='b';
        case 'avg_abg'  
            GA_str='c';
        case 'avg_abg_pa'
            GA_str='d';
            
        case '2a-b_pa'    
            GA_str='2';
            
        case '2ab_a_b'
            GA_str='e';
            
        case 'b_pa'    
            GA_str='bP';
            
        case 'none'
            GA_str='NA';
            
        otherwise
            GA_str='Z';
    end
    
    if strcmp(params.breed_mode1,'none') == 0

        switch params.GA_metric

            case 'chi'
                GA_str=[GA_str,'C'];

            case 'sharpness'
                GA_str=[GA_str,'S'];

            case 'area'
                GA_str=[GA_str,'A'];

            case 'summed_phase'
                GA_str=[GA_str,'P'];
                
            case 'TV'
                GA_str=[GA_str,'V'];

        end
        
    end
    
    npop=num2str(params.population);
    ngen=num2str(params.generations);
    
    GA_return=upper(params.GA_return);
    
    %if strcmp(GA_return,'BEST'), GA_return='BST';end
    switch GA_return
        
        case {'BEST','BST'}
            GA_return='BST';
        
        case {'BEST-CULL'}
            GA_return='BSTc';
            
        
        case {'AVG-HALF-CULL'}
            GA_return='AVGhc';    
            
        case {'AVG','AVERAGE'}
            GA_return='AVG';
            
        case {'AVG-HALF','AVERAGE-HALF'}
            GA_return='AVGh';
            
    end
    
    name1=['GA',ALG2,'-',GA_return,'-',GA_str,'-',npop,'-',ngen,'-',its,'-',pc,'-',son];

else
    name1=[ALG1,ALG2,its,'-',pc,'-',son];
end
name=[name,name1];
    
    
if params.do_2D == 1, name=['2D-',name];end

if params.iterate_avg_enabled == 1,name=[name,'-AVG'];end

end

