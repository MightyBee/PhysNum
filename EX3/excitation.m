% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'excitation';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

theta0=0;
thetadot0=1e-2;
Omega=2*sqrt(9.81/0.1);
d=0.005;
kappa=0.05;
dt=0.002;
tFin=100;

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = [repertoireOut,filename,'_kappa=',num2str(kappa),'_dt=',num2str(dt), '.out'];
cmd = sprintf('%s%s %s tFin=%.15g d=%.15g Omega=%.15g theta0=%.15g thetadot0=%.15g dt=%.15g output=%s', repertoireExe, executable, input, tFin, d, Omega, theta0, thetadot0, dt, outputFile);
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

W=zeros(size(P,1),1);
for i = 2:size(P,1)
    W(i) = W(i-1) + (P(i-1)+P(i))*dt/2;
end


% Figures
figure('Position',[50,50,600,400]);
plot(t,energy,t, W,'k.',t,energy-W, 'r.')
xlabel('t [s]')
ylabel('Energy [J]')
set(gca,'fontsize',15);
grid on
print('figures/energy&W', '-depsc');

