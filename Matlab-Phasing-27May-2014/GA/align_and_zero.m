function pn=align_and_zero(ref,pn)

%check the conj reflection
pn=check_conj_ref(ref,pn);

%align just the amplitudes intiailly
%shift iterates to the dominant one
[h k l]=register_3d_reconstruction(abs(ref),abs(pn));
pn=zero_phase(((sub_pixel_shift((pn),h,k,l))));
%pn=zero_phase(pn);


end