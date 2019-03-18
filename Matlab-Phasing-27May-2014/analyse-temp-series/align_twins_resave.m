function [ output_args ] = align_twins_resave(top_dir,prefix)
%jclark
% checks the saved reconstructions 
% and resaves if the twin is present.
% used for temperature series
% eg. align_twins_resave(dir,'Rec-6x8-*-GAHIO-AVGh-mC-10-5-200-')
% will align all the reconstructions with the above prefix, the * allows
% different file numbers.
% will only resave if it is the twin
% if you want all to be resaved regardless use align_vtk_resave.m


alpha=1;        %which reconstruction to align to

try
    top_dir;
catch
    name_of_this_file='analyse_twins_resave';
    dir_file=which(name_of_this_file);    %finds the current folder for the phasing
    top_dir=dir_file(1:findstr(dir_file,name_of_this_file)-1);   %keep just the directory
end

try
    prefix;
catch
    prefix='rand-starts/CVL/CVL';
end

%%
% get names of the summed and averaged reconstructions and a paramter file
% for each reconstruction


%string='rand-starts/CVL/CVL-AMP.rec';
string=[prefix,'**-AMP.rec'];
amp_fs=rdir([top_dir,'**/*/',string]);

%string='rand-starts/CVL/CVL-PH.rec';
string=[prefix,'**-PH.rec'];
ph_fs=rdir([top_dir,'**/*/',string]);

%get the lab frame matlab rec if it exists
string=[prefix,'**-LAB.rec'];
lab_fs=rdir([top_dir,'**/*/',string]);
string=[prefix,'**-LAB-S.rec'];
labs_fs=rdir([top_dir,'**/*/',string]);
if size(amp_fs,1) == size(lab_fs,1),output_lab=1;else output_lab=0;end


n_temps=size(amp_fs,1);


first = load_rec(strtrim(char(amp_fs(alpha).name)),strtrim(char(ph_fs(alpha).name)));


for qq = 1:n_temps
    
    reflected=0;
    
    if qq ~= alpha
        
        pn = load_rec(strtrim(char(amp_fs(qq).name)),strtrim(char(ph_fs(qq).name)));

        try
            temp=pn;
            if sum(size(temp)-size(first)) ~= 0
                disp('Resizing for comparison....')
                tempFFt=fftshift(fftn(fftshift(temp)));
                tempFFt=zero_pad_ver3(tempFFt,size(first,2),size(first,1),size(first,3));
                temp=fftshift(ifftn(fftshift(tempFFt)));
                
            end
            
            cnj_rf=is_conj_ref(abs(first),abs(temp));

            if cnj_rf ~=0,
                second=conj_reflect(pn);
                disp([num2str(qq), ' - REFLECTED!!!'])
                reflected=1;
            else
                second=(pn);
            end

            array=abs(second);
            save(amp_fs(qq).name,'array')

            array=atan2( imag(second),real(second) );
            save(ph_fs(qq).name,'array')
            
            if reflected == 1 && output_lab == 1
               save_dir=extract_dir_from_string( strtrim(char(lab_fs(qq).name)) );
               disp('Resaving lab output and .vtk....') 
               load(strtrim(char(lab_fs(qq).name)),'-mat')
               array=conj_reflect(array);
               save(strtrim(char(lab_fs(qq).name)),'array')
               savevtk2scalar(abs(array),[save_dir,'Amp-Phase.vtk'],angle(array),1)
               
               load(strtrim(char(labs_fs(qq).name)),'-mat')
               array=conj_reflect(array);
               save(strtrim(char(labs_fs(qq).name)),'array')
               
                   
               savevtk2scalar(abs(array),[save_dir,'Support.vtk'],[],1)  
            end
            
        catch
            disp('Error with ....')
            disp(strtrim(char(amp_fs(qq).name)))
        end
        pn=0;
        array=0;
        
        second=0;
   end
    
    
    
end


end

function cnj_rf=is_conj_ref(a,b)

c1=cross_correlation(a,conj_reflect(b));
c2=cross_correlation(a,b);

%c1=cross_correlation(a,b);
%c2=cross_correlation(a,b);


c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end
