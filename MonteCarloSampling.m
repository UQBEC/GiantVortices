clear
close all
clc

%% Ellipse and Vortex Parameters
xi = 0.53;
a = 120/2/xi; 
b = 85/2/xi; 
% m the parameter, see http://mathworld.wolfram.com/Parameter.html
% k = sqrt(m) is the Elliptic Modulus, see http://mathworld.wolfram.com/EllipticModulus.html
q = ((a - b)/(a + b))^2;
k = 1;
for n = 1:100
    k = k*((1 + q^(2*n))/(1 + q^(2*n-1)) )^8;
end
k = sqrt(16*q*k);
m = k^2;

N = 18;
kappa = [ones(N/2,1); -ones(N/2,1)];

%% Check the conformal map by applying to some concentric ellipses
phi = 0:pi/100:2*pi;
for lambda = 0.1:0.1:1
    temp = lambda*(a*cos(phi) + b*1i*sin(phi));
    figure(1)
    subplot(121)
    plot(temp)
    axis square
    xlim([-a a]*1.2)
    ylim(xlim)
    hold on
    [temp,~] = ellipse2circlemap(temp,m,a,b);
    subplot(122)
    plot(temp)
    axis square
    xlim([-1.2 1.2])
    ylim(xlim)
    hold on
end
subplot(121)
xlabel('$x$','interpreter','Latex','Fontsize',18)
ylabel('$y$','interpreter','Latex','Fontsize',18)
title('$z \in \Omega$','interpreter','Latex','Fontsize',18)
set(gca,'Fontsize',14)
subplot(122)
xlabel('$\zeta_x$','interpreter','Latex','Fontsize',18)
ylabel('$\zeta_y$','interpreter','Latex','Fontsize',18)
title('$\zeta \in \mathcal{D}$','interpreter','Latex','Fontsize',18)
set(gca,'Fontsize',14)
drawnow

%% Do 10 runs of Monte Carlo with 1e8 Samples
%Nsamples = 1e8;

Nsamples = 1e6;
for runs = 1%:10
H = zeros(Nsamples,1);
D = zeros(Nsamples,1);
ell = zeros(Nsamples,1);
%parpool(4)
tic
parfor (zz = 1:Nsamples,4)
 if floor (zz/10000) == zz/10000
     disp(num2str(zz))
 end

%Generate positions
r = sqrt(rand(N,1));
theta = 2*pi*rand(N,1);
x = a*r.*cos(theta);
y = b*r.*sin(theta);
z = x+1i*y;

%Calculate energy (H), complex dipole moment (D = Dx + iDy) and n.n
%distance (ell). Note Dx is the ordrer parameter for clustering in the
%ellipse
[zeta,dzeta] = ellipse2circlemap(z,m,a,b);
H(zz) = ConformalEnergy(zeta,dzeta,kappa);
D(zz) = kappa.'*z;
ell(zz) = getNearestNeighbourDistance(x,y,N);

% Unfold this cell to plot the vortex positions
%{
figure(2)
clf
subplot(121)
plot(z(1:end/2),'.r')
hold on
plot(z(end/2+1:end),'.b')
plot(a*cos(phi) + 1i*b*sin(phi),'-k')
axis square
xlim([-150 150])
ylim(xlim)
subplot(122)
plot(zeta(1:end/2),'.r')
hold on
plot(zeta(end/2+1:end),'.b')
plot(exp(1i*phi),'-k')
axis square
xlim([-1.2 1.2])
ylim(xlim)
pause(0.1)
%}
end
toc

%Calculate density of states and normalize (here normalized to unity) 
[counts,energies] = hist(H/N,500);
dE = energies(2)- energies(1);
[~,ind] = max(counts);
Em = energies(ind);
counts = counts/Nsamples/dE;
figure(1)
plot(energies-Em,counts,'k')

save(['a120_b85_N18_Nsamples' num2str(Nsamples) '_D_ell_run' num2str(runs) '.mat'],'H','energies','counts','D','ell','N')

end
