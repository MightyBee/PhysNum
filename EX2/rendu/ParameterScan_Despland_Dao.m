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

nsimul = 30; % Nombre de simulations a faire

nsteps = round(logspace(2.1,4.95,nsimul)); % Nombre d'iterations entier de 10^2 a 10^4
tfin = 1.09e-7; % TODO: Remplacer la valeur de tfin
dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner
param = nsteps; % Valeurs du parametre a scanner



%% Simulations %%
%%%%%%%%%%%%%%%%%
schemas={'E','EC','RK2','V'};
nschemas=size(schemas,2);
%cmd = sprintf('rm -r simulations && mkdir simulations && mkdir simulations/%s && mkdir simulations/%s && mkdir simulations/%s',schemas{1},schemas{2},schemas{3})
%disp(cmd)
%system(cmd);
output = cell(nschemas, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for k = 1:nschemas
  for i = 1:nsimul
    output{k,i} = ['simulations/',schemas{k},'/',paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s schema=%s %s=%.15g output=%s', repertoire, executable, input, schemas{k}, paramstr, param(i), output{k,i});
    disp(cmd)
    system(cmd);
  end
end
%% Analyse %%
%%%%%%%%%%%%%
%filename = 'configuration.in';
%delimiterIn = '=';
%headerlinesIn = 0;
%A = importdata(filename,delimiterIn,headerlinesIn);
%disp(A);
%disp(A(1,1).data(5,1));
v0=4e5;
omega=(3*1.6022/1.6726)*1e8;
x0=-1.39e-3;%-v0/omega;

error = zeros(nschemas,nsimul);
for k = 1:nschemas
  for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{k,i}); % Chargement du fichier de sortie de la i-ieme simulation
    t = data(:,1);
    x = data(:,2);
    y = data(:,3);
    x_th = v0/omega*(1-cos(omega*t))+x0; % TODO: Entrer la vraie solution analytique en fonction du temps
    y_th = v0/omega*sin(omega*t); % TODO: Entrer la vraie solution analytique en fonction du temps
    error(k,i) = max(sqrt((x-x_th).^2+(y-y_th).^2));
  end
end

figure
loglog(dt, error(1,:),'r+', dt, error(2,:),'b+', dt, error(3,:), 'k+', dt , error(4,:), 'm+')
xlabel(['\Delta','t [s]'], 'fontsize', 15)
ylabel('Maximum de l''erreur sur la position [m]', 'fontsize', 15)
set(gca,'fontsize',15)
lgd = legend('Euler', 'EulerCromer', 'RungeKutta2','Verlet', 'Location', 'southeast');
set(lgd,'FontSize', 15);
grid on


print('etudeConv', '-depsc');
%set(gcf, 'Units', 'inches');
%screenposition = get(gcf, 'Position');
%set(gcf, 'PaperPosition', [0 0 screenposition(3:4)],
%         'PaperSize', [screenposition(3:4)]);
%print('etudeConv', '-dpdf');
