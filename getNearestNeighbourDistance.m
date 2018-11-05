function ell = getNearestNeighbourDistance(x,y,N)
Rij = zeros(N); %Matrix of distances
for ii = 1:N
    Rij(:,ii) = sqrt( (x-x(ii)).^2 + (y-y(ii)).^2 );
    Rij(ii,ii) = NaN; %exclude diagonal terms
end
ell = mean(nanmin(Rij));
return

