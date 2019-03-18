function align_vtk_resave( top_dir,prefix )
%jclark

conj_first=1;
flipdim=1;

%get the lab frame matlab rec if it exists
string=[prefix,'**-LAB.rec'];
lab_fs=rdir([top_dir,'**/',string]);
string=[prefix,'**-LAB-S.rec'];
labs_fs=rdir([top_dir,'**/',string]);

output_lab=1;
n_temps=size(lab_fs,1);

%load first rec
load(strtrim(char(lab_fs(1).name)),'-mat');
first=array;
array=[];
%get its dir
first_dir=extract_dir_from_string(strtrim(char(lab_fs(1).name)));
%get spacing from existing vtk
spacing=get_spacing_vtk([first_dir,'Amp-Phase.vtk']);

%do we want to relfect the first?
if conj_first == 1,
    disp('Conjugating and reflecting first....')
    first=conj_reflect(first);end

%zero the phase
first=zero_phase(first,0);

%flip the dims for the vtk output
if flipdim == 1,first=flip_all_dim(first);end

%resave the first one
savevtk2scalar(abs(first),[first_dir,'Amp-Phase.vtk'],angle(first),spacing);

if n_temps ~= 1
    for qq = 2:n_temps

        reflected=0;

        %load the others
        load(strtrim(char(lab_fs(qq).name)),'-mat');
        pn=array;
        array=[];
        temp=pn;

        save_dir=extract_dir_from_string(strtrim(char(lab_fs(qq).name)));

        %check if it is the same size, prob not necessary at this stage
        if sum(size(temp)-size(first)) ~= 0
            disp('Resizing for comparison....')
            tempFFt=fftshift(fftn(fftshift(temp)));
            tempFFt=zero_pad_ver3(tempFFt,size(first,2),size(first,1),size(first,3));
            temp=fftshift(ifftn(fftshift(tempFFt)));
        end

        %is the curent one reflected relative to the first?
        cnj_rf=is_conj_ref(abs(first),abs(temp));

        if cnj_rf ~=0,
            second=conj_reflect(pn);
            disp([num2str(qq), ' - REFLECTED!!!'])
            reflected=1;
        else
            second=(pn);
        end
        pn=[];
        second=zero_phase(second);

        if flipdim == 1,second=flip_all_dim(second);end
        savevtk2scalar(abs(second),[save_dir,'Amp-Phase.vtk'],angle(second),spacing)

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
