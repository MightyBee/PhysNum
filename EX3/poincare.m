% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'poincare';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

sampling=100;
theta0=0;
thetadot0=1e-1;
Omega=15.1198%sqrt(9.81/0.1); %%% 15.1197 pour theta0=0, thetadot0=1e-1 et d=0.04
d=0.04;
kappa=0;
dt=2*pi/(sampling*Omega);
tFin=dt*10000*sampling;

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = [repertoireOut,filename,'_theta0=',num2str(theta0),'_thetadot0=',num2str(thetadot0), '.out'];
cmd = sprintf('%s%s %s tFin=%.15g d=%.15g Omega=%.15g theta0=%.15g thetadot0=%.15g dt=%.15g output=%s sampling=%d', repertoireExe, executable, input, tFin, d, Omega, theta0, thetadot0, dt, outputFile, sampling);
disp(cmd)
system(cmd);

% Chargement des donnees
output = load(outputFile);

% Extraction des quantites d'interet
t = output(:,1);
theta = output(:,2);
thetadot = output(:,3);
energy= output(:,4);
P = output(:,5);
clear output

% Figures
figure('Position',[50,50,600,400]);
plot(theta, thetadot, 'k.')
xlabel('\theta [rad]')
ylabel('\omega [rad/s]')
set(gca,'fontsize',15);
grid on
print(['figures/poincare'], '-depsc');
