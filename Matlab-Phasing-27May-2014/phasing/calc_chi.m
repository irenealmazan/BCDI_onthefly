function error = calc_chi(Mk,Mm)
% calculates error
nume=sum(sum(sum(abs(Mk-Mm).^2)));                  %calculating the error
denom=sum(sum(sum(abs(Mm).^2)));
error=nume/denom;

end

