% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'excitation';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

theta0=0;
thetadot0=1e-2;
w0=sqrt(9.81/0.1);
Omega=1*w0;
d=0.03;
kappa=0;
dt=0.002;
tFin=250;

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = [repertoireOut,filename,'_kappa=',num2str(kappa),'_dt=',num2str(dt), '.out'];
cmd = sprintf('%s%s %s tFin=%.15g d=%.15g Omega=%.15g kappa=%.15g theta0=%.15g thetadot0=%.15g dt=%.15g output=%s', repertoireExe, executable, input, tFin, d, Omega, kappa, theta0, thetadot0, dt, outputFile);
disp(cmd)
system(cmd);

% Chargement des donnees
output = load(outputFile);

% Extraction des quantites d'interet
t = output(:,1);
energy= output(:,4);
P = output(:,5);
clear output

W=zeros(size(P,1),1);
for i = 2:size(P,1)
    W(i) = W(i-1) + (P(i-1)+P(i))*dt/2;
end


% Figures
fig1=figure('Position',[50,50,600,400]);
plot(t,energy,t, W,'k.',t,energy-W, 'r.')
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
xlabel('t [s]')
ylabel('Energy [J]')
xlim([0 tFin])
ylim([-0.005 0.04]);
set(gca,'fontsize',15);
lgd = legend('$E_{mec}$', '$W_{\rm{NC}}$','$E_{mec}-W_{\rm{NC}}$', 'Location', 'northwest','Interpreter','latex');
set(lgd,'FontSize', 13);
grid on
print(fig1,'figures/energyW', '-depsc');

indexOfInterest = (t < 95) & (t > 90); % range of t near perturbation
fig2=figure('Position',[50,50,600,400]);
plot(t(indexOfInterest),energy(indexOfInterest),t(indexOfInterest),W(indexOfInterest),'k.','Linewidth', 1.5) % plot on new axes
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
xlabel('t [s]')
ylabel('Energy [J]')
set(gca,'fontsize',15);
lgd = legend('$E_{mec}$', '$W_{\rm{NC}}$', 'Location', 'northeast','Interpreter','latex');
set(lgd,'FontSize', 13);
print(fig2,'figures/energyWZoom', '-depsc');

fig3=figure('Position',[50,50,600,400]);
plot(t,energy-W, 'r.')
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
xlabel('t [s]')
ylabel('Energy [J]')
xlim([0 tFin])
set(gca,'fontsize',15);
lgd = legend('$E_{mec}-W_{\rm{NC}}$', 'Location', 'northwest','Interpreter','latex');
set(lgd,'FontSize', 13);
grid on
print(fig1,'figures/energyW', '-depsc');

