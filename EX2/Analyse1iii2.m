% Nom du fichier d'output a analyser
repertoireOut = 'stabilite';
schema = {'Euler','EulerCromer','RungeKutta2'};
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base


nsteps = 500; % Nombre d'iterations entier de 10^2 a 10^4
nsimul = 3; % Nombre de simulations a faire
tfin = 1.09e-7; % TODO: Remplacer la valeur de tfin
dt = tfin ./ nsteps;

paramstr = 'schema'; % Nom du parametre a scanner
param = schema; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
  output{i} = [repertoireOut,'/',paramstr, '=', param{i}, '.out'];
  % Execution du programme en lui envoyant la valeur a scanner en argument
  cmd = sprintf('%s%s %s nsteps=80 %s=%s output=%s', repertoireExe, executable, input, paramstr, param{i}, output{i});
  disp(cmd)
  system(cmd);
end


output1 = load(output{1});
output2 = load(output{2});
output3 = load(output{3});

% Extraction des quantites d'interet
t1 = output1(:,1);
x1 = output1(:,2);
y1 = output1(:,3);
vx1 = output1(:,4);
vy1 = output1(:,5);
energy1 = output1(:,6);
mu1 = output1(:,7);

t2 = output2(:,1);
x2 = output2(:,2);
y2 = output2(:,3);
vx2 = output2(:,4);
vy2 = output2(:,5);
energy2 = output2(:,6);
mu2 = output2(:,7);

t3 = output3(:,1);
x3 = output3(:,2);
y3 = output3(:,3);
vx3 = output3(:,4);
vy3 = output3(:,5);
energy3 = output3(:,6);

clear output1
clear output2
clear output3


fig1=figure('Position',[20,20,700,400])
semilogy(t1,energy1,t2,energy2,t3,energy3)
grid on
xlabel('t [s]','Fontsize', 13)
ylabel('E [J]','Fontsize', 13)
lgd4=legend('Euler', 'Euler-Comer','Runge-Kutta 2','Location','northwest')
set(lgd4,'FontSize', 12);
set(gca,'fontsize',13);

print(fig1,'etudeStabiliteCompEnergy', '-depsc');

%print([repertoireOut,'/',filename], '-depsc');
