function pn=check_conj_ref(ref,pn)

cnj_rf=is_conj_ref(abs(ref),abs(pn));
            
if cnj_rf ~=0,pn=conj_reflect(pn);end

end
