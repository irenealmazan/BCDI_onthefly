function [ hparams ] = calc_cc_shells(hparams)
%jclark

nthings=max(size(hparams.amp_shell_vals));  %get number of objects

first_amp=hparams.amp_shell_vals{1};

first_ph=hparams.ph_shell_vals{1};

nshells=max(size(first_amp));

hparams.nthings=nthings;
hparams.nshells=nshells;

for qq = 1:nthings
    
    second_amp=hparams.amp_shell_vals{qq}; %get second shell values
    second_ph=hparams.ph_shell_vals{qq};
    
    for ww = 1:nshells
        
        one_amp=first_amp{ww}; %get shell values
        one_ph=first_ph{ww};

        two_amp=second_amp{ww}; %get shell vals from second
        two_ph=second_ph{ww}; %get shell vals from second
        
        
        hparams.cc_amp(qq,ww)=max(max(normxcorr2(one_amp,two_amp)));
        hparams.cc_ph(qq,ww)=max(max(normxcorr2(one_ph,two_ph)));
       
    end
    
    
end

end

