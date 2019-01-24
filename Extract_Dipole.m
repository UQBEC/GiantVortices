clear all
close all
clc
%%
% M. T. Reeves 01/11/2018

% This script generates a 2D energy-dipole moment histogram, and extracts
% the most likely value of the dipole moment from the peak of the
% distribution at each energy.

% !Warning: you need a lot of RAM to run this as there are 10 samples each
% containing 10^8 values. 

Em = 3.1547;  %Energy at peak entropy
N = 18;       %Vortex Number
a = 60/0.53;  %Ellipse major axis / healing length (in um)
Ds = 0.47*a;  %upper limit of dipole moment (supercondensation)
Nruns = 10;

%Energy and dipole moment vectors for 2D histogram
Dvals = linspace(-1.4,1.4,60);
Evals = linspace(-2.5,8.4,50);
dE = Evals(2) - Evals(1);
dD = Dvals(2) - Dvals(1);
ddd = [];

% Load each run of Monte Carlo data and generate E-D histogram. Extract
% peak at each value of E to determine the most likely value of D. 
for runs = 1:Nruns
    disp(num2str(runs))
    load(['a120_b85_N18_Nsamples1e8_.100pc_D_ell_CB_rad_run' num2str(runs) '_.mat']);
    %load(['a120_b85_N18_Nsamples1e6_D_ell_run' num2str(runs) '.mat']);
    figure(runs)
    EDhist =  histogram2(H/N-Em,abs(real(D))/Ds/N,Evals,Dvals);
    eval(['test' num2str(runs) ' = EDhist;'])
    [derp,ind] = max(EDhist.Values,[],2);
    dd = EDhist.YBinEdges(ind);
    ddd = [ddd; dd]; %Collect values for each run.
end

%% Plotting
%Combine the different runs into single histogram

%Plot
figure(666)
clf
hold on
plot(Evals(1:end-1)-dE/2,mean(abs(ddd),1)-dD/2,'ok','Linewidth',1.0,'MarkerFaceColor',[0.7 0.7 0.7])
axis square
xlim([-1.5 3])
xlabel('Energy per vortex, $(E-E_m)/N$','Interpreter','Latex','Fontsize',18)
ylabel('Dipole moment, $D/D_s$','Interpreter','Latex','Fontsize',18)
%% Fit to compare with mean-field prediction
xdata = (Evals(1:end-1)-dE/2).';
ydata = (mean(abs(ddd),1)-dD/2).';

%Ignore noisy data at high and low energies where sample numbers are low
ind = (xdata > -1.5 ) & (xdata < 3);
xdata = xdata(ind);
ydata = ydata(ind);
fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0,0],...
               'Upper',[1,1],...
               'StartPoint',[1 1]);
%Mean-field curve is 0 below Ec and sqrt(E) above Ec. y = real(sqrt(x-a)) 
%conveniently gives equivalent to piecewise fit.

%note that values obtained for D0 and Ec exhibit some dependence on the
%choice of binning (a few % ).
ft = fittype('D0*real((x-Ec).^0.5)','options',fo);
[curve2,gof2] = fit(xdata,ydata,ft);
figure(666)
hold on
xdata2 = linspace(-1.5,3,200)
plot(xdata2,curve2(xdata2),'-','Color',[0.5 0 0.5],'Linewidth',1)
hold on
plot(curve2.Ec,0,'p','Color',[0.5 0 0.5],'MarkerFaceColor',[0.85 0.7 1],'Markersize',14)









