function H = ConformalEnergy(zeta,dzeta,kappa)

%Calculates the energy under conformal mapping to the unit circle.
    
N = length(kappa);
Hself = -sum(kappa.^2.*log(abs(dzeta./(1-abs(zeta).^2))));

Hint = 0;
for jj = 1:N
    fij = abs( (zeta(jj) - zeta)./(1 - zeta(jj)*conj(zeta)) );
    fij = kappa(jj)*kappa.*log(fij);
    fij(jj) = 0;
    Hint = Hint - sum(fij);
end
H = (Hself + Hint);

end


