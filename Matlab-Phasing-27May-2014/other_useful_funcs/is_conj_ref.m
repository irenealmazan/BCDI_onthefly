function cnj_rf=is_conj_ref(a,b)

c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));

c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end