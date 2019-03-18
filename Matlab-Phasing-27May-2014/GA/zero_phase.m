function pn=zero_phase(pn,val)

    try
        val;
    catch
        val=0;
    end

    ph=atan2(imag(pn),real(pn));

    SS=shrink_wrap(abs(pn),.2,.5);  %get just the crystal, i.e very tight support
    avg_ph=mean(ph(find(SS > 0)));
    ph=ph-avg_ph+val;

    pn=abs(pn).*exp(i*ph);


end