clear 
clc
close all

Nruns = 10;
xi = 0.53;      % Average healing length
xi_sim = 0.5;   % Healing length at peak density

a = 60/xi_sim;   %Major axis
b = 42.5/xi_sim; %Minor axis
Ds = 0.47*a;     %Dipole moment at supercondensation
DT = 20;         %time spacing between simulation outputs (in sim units)
delta = 5;       %data only extracted for every 5th output
t0 = 0.5/(1460); %time unit for simulations
Em = 3.154;      %Energy at peak entropy

%% Paddle
loaddir = './PaddleRuns/';      %load directory
load([loaddir 'run1.mat'])      %load the data
t = (1:length(H))*DT*t0*delta;  % time vector for plotting

[~,ind] = (min(abs(t-1))); %find index where t = 1s
H_plot = H;
ell_plot = ell_sim_comb;
N_plot = N_sim_comb;      
Dx_plot = (Dx_paddle);
Dx2_plot = (Dx2_paddle);
ell0 = sqrt(a*b./N_sim_comb)*sqrt(0.89);
for run = 2:10
    load([loaddir 'run' num2str(run) '.mat'])      
    H_plot = H_plot + H;
    ell_plot = ell_plot + ell_sim_comb;
    N_plot = N_plot + N_sim_comb;      
    Dx_plot = Dx_plot + (Dx_paddle);
    Dx2_plot = Dx2_plot + (Dx2_paddle);
    ell0 = ell0 + sqrt(a*b./N_sim_comb)*sqrt(0.89);
    
end

%rolling average over 1s -- (note time average is required for fair
%comparison with microcanonical ensemble, to average over fluctuations/oscillations
%in dipole moment under dynamical evolution. (D is not conserved under
%dynamical evolution).
H_plot = movmean(H_plot,ind)/Nruns;
ell_plot = movmean(ell_plot./ell0,ind)/Nruns;
N_plot = movmean(N_plot,ind)/Nruns;
Dx_plot = (movmean(abs(Dx_plot),ind))/Nruns;
Dx2_plot = (movmean(abs(Dx2_plot),ind))/Nruns;


figure()
subplot(221)
plot(t,N_sim_comb,'Linewidth',2)
ylim([0 40])
xlim([0 t(end)])
xlabel('Time [s]','Interpreter','Latex')
ylabel('Vortex Number, $N$','Interpreter','Latex')
subplot(222)
plot(t,H_plot./N_plot-Em,'Linewidth',2)
ylabel('Energy, $(E-E_m)/E_0 N$','Interpreter','Latex')
xlabel('Time [s]','Interpreter','Latex')
xlim([0 t(end)])
subplot(223)
plot(t,ell_plot,'Linewidth',2)
xlabel('Time [s]','Interpreter','Latex')
ylabel('mean nearest-neighbour distance, $\ell/\ell_0$','Interpreter','Latex')
xlim([0 t(end)])
subplot(224)
plot(t,abs(Dx_plot)/Ds./N_plot(1),'Linewidth',2)
hold on 
plot(t,abs(Dx2_plot)/Ds./N_plot(1),'--','Linewidth',2,'Color',[0.7 0.7 1])
xlabel('Time [s]','Interpreter','Latex')
ylabel('Dipole moment, $D/D_s$','Interpreter','Latex')
ylim([0 1])
xlim([0 t(end)])


%% Comb

clear 
clc

Nruns = 10;
xi = 0.53;
xi_sim = 0.5;

a = 60/xi_sim;   %Major axis
b = 42.5/xi_sim; %Minor axis
Ds = 0.47*a;     %Dipole moment at supercondensation
DT = 20;         %time spacing between simulation outputs (in sim units)
delta = 5;       %data only extracted for every 5th output
t0 = 0.5/(1460); %time unit for simulations
Em = 3.154;      %Energy at peak entropy

loaddir = './CombRuns/';      %load directory
load([loaddir 'run1.mat'])      %load the data
t = (1:length(H))*DT*t0*delta;  % time vector for plotting

[~,ind] = (min(abs(t-1))); %find index where t = 1s
H_plot = H;
ell_plot = ell_sim_comb;
N_plot = N_sim_comb;      
Dx_plot = abs(Dx_comb);
ell0 = sqrt(a*b./N_sim_comb)*sqrt(0.89);
for run = 2:10
    load([loaddir 'run' num2str(run) '.mat'])      
    H_plot = H_plot + H;
    ell_plot = ell_plot + ell_sim_comb;
    N_plot = N_plot + N_sim_comb;      
    Dx_plot = Dx_plot + abs(Dx_comb);
    ell0 = ell0 + sqrt(a*b./N_sim_comb)*sqrt(0.89);
    
end
H_plot = movmean(H_plot,ind)/Nruns; 
Dx_plot = movmean(abs(Dx_plot),ind)/Nruns;
N_plot = movmean(N_plot,ind)/Nruns;
ell_plot = smooth(ell_plot./ell0,ind);

    
figure(gcf)
subplot(221)
hold on
plot(t,N_plot,'Linewidth',2,'Color',[1 0.5 0])
subplot(222)
hold on
plot(t,H_plot./N_plot-Em,'Linewidth',2,'Color',[1 0.5 0])
hold on
plot(xlim,[0 0],'-.k','Linewidth',1)
plot(xlim,[1 1]*0.81,'--','Linewidth',1,'Color',[0.5 0 1])
xlim([0 t(end)])
subplot(223)
hold on
plot(t,ell_plot,'Linewidth',2,'Color',[1 0.5 0])
plot(xlim,[1 1],'--k')
xlim([0 t(end)])
subplot(224)
hold on
plot(t,(Dx_plot)/Ds./N_plot(1),'Linewidth',2,'Color',[1 0.5 0])
plot(xlim,[1 1],'--r','Linewidth',1)
ylim([0 1.2])
xlim([0 t(end)])


%% Paddle histograms

loaddir = './PaddleRuns/';                %load directory
load([loaddir 'run1Histograms.mat'])      %load the data
OMEGAP = omegap;
OMEGAM = omegam;
for run = 2:10
    load([loaddir 'run' num2str(run) 'Histograms.mat'])      %load the data
    OMEGAP = OMEGAP + omegap;
    OMEGAM = OMEGAM + omegam;
end

figure
OMEGA_MAX = max(abs(OMEGAP(:) - OMEGAM(:)));
OMEGA = (OMEGAP-OMEGAM)./OMEGA_MAX
imagesc(OMEGA)

color1 = [1 0 0];
color2= [0 0 1];
color3 = [1 1 1];

map = [linspace(color1(1),color3(1),128).' linspace(color1(2),color3(2),128).' linspace(color1(3),color3(3),128).']
map = [map; linspace(color3(1),color2(1),128).' linspace(color3(2),color2(2),128).' linspace(color3(3),color2(3),128).']
caxis([-1 1])
colormap(map)

%%
loaddir = './CombRuns/';                %load directory
load([loaddir 'run1Histograms.mat'])      %load the data
OMEGAP = omegap;
OMEGAM = omegam;
for run = 2:10
    load([loaddir 'run' num2str(run) 'Histograms.mat'])      %load the data
    OMEGAP = OMEGAP + omegap;
    OMEGAM = OMEGAM + omegam;
end

figure
OMEGA = (OMEGAP-OMEGAM);
OMEGA = OMEGA./OMEGA_MAX;
imagesc(OMEGA)

color1 = [1 0 0];
color2= [0 0 1];
color3 = [1 1 1];

map = [linspace(color1(1),color3(1),128).' linspace(color1(2),color3(2),128).' linspace(color1(3),color3(3),128).']
map = [map; linspace(color3(1),color2(1),128).' linspace(color3(2),color2(2),128).' linspace(color3(3),color2(3),128).']
caxis([-1 1])
colormap(map)








    
    