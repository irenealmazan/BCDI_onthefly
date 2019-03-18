function constrained = phase_gradient_constraint(pn,phx )
%jclark

dims=ndims(pn);

expphi=exp(i*angle(pn));
%expphi=pn./abs(pn);

if dims == 2
    
    [ fx fy] = phase_gradient( expphi,[2,2,2]);

    fx(abs(fx) <= phx) = phx; %constrain the gradient
    
    const=zeros(size(pn))+phx;  %create a constant array
    
    %phase=angle(pn).*(const./abs(fx));%abs(fx)./const;
    
    %constrained=abs(pn).*exp(i*phase);
    
    constrained=abs(pn).*(expphi).^(const./abs(fx));
    %constrained=abs(pn).*expphi.*exp(angle(pn).*((const./abs(fx)-1)));
    
end





end

