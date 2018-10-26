% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 10; % Nombre de simulations a faire

nsteps = round(logspace(2.2,3.5,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4
tfin = 1.09e-7; % TODO: Remplacer la valeur de tfin
dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner
param = nsteps; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
  output{i} = ['application2/simulations/',paramstr, '=', num2str(param(i)), '.out'];
  % Execution du programme en lui envoyant la valeur a scanner en argument
  cmd = sprintf('%s%s %s Kappa=100 %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
  disp(cmd)
  system(cmd);
end
%% Analyse %%
%%%%%%%%%%%%%

pFinal = zeros(2,nsimul);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
  data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
  pFinal(1,i) = data(end,2);
  pFinal(2,i) = data(end,3);
end

dtSquared= dt.*dt;
fit1= polyfit(dtSquared, pFinal(1,:),1);
fit2= polyfit(dtSquared, pFinal(2,:),1);
dtFit= linspace(0,(tfin/nsteps(1))^2);
x1 = polyval(fit1,dtFit);
x2 = polyval(fit2,dtFit);

figure('Position',[20,20,650,400])
plot(dtSquared, pFinal(1,:),'r+', dtSquared, pFinal(2,:),'b+');
hold on
plot(dtFit,x1,'--r',dtFit,x2,'--b')
xlabel(['(\Delta','t)^2 [s^2]'], 'fontsize', 13)
ylabel('Position finale [m]', 'fontsize', 13)
set(gca,'fontsize',13)
lgd = legend('x_{final}', 'y_{final}', 'Location', 'southeast');
set(lgd,'FontSize', 13);
%hold off
grid on

print('application2/etudeConv2', '-depsc');
%set(gcf, 'Units', 'inches');
%screenposition = get(gcf, 'Position');
%set(gcf, 'PaperPosition', [0 0 screenposition(3:4)],
%         'PaperSize', [screenposition(3:4)]);
%print('application1/etudeConv1', '-dpdf');

figure
subplot(1,2,1)
semilogx(dt, pFinal(1,:),'r+');
xlabel(['\Delta','t [s]'], 'fontsize', 13)
ylabel('x_{final} [m]', 'fontsize', 13)
set(gca,'fontsize',13)
grid on

subplot(1,2,2)
semilogx(dt, pFinal(2,:),'b+');
xlabel(['\Delta','t [s]'], 'fontsize', 13)
ylabel('y_{final} [m]', 'fontsize', 13)
set(gca,'fontsize',13)
grid on


print('application2/etudeConv6', '-depsc');
