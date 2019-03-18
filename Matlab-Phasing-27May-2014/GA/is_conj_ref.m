function cnj_rf=is_conj_ref(a,b,nos)

try
    nos;
catch
    nos=0;
end

if nos == 0
    c1=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(conj_reflect(b),.1,.1));
    c2=cross_correlation(shrink_wrap(a,.1,.1),shrink_wrap(b,.1,.1));
else
    c1=cross_correlation(a,conj_reflect(b));
    c2=cross_correlation(a,b);
end
%c1=cross_correlation(a,b);
%c2=cross_correlation(a,b);


c1=max(c1(:));
c2=max(c2(:));

if c1 > c2,cnj_rf=1;else cnj_rf=0;end


end