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
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 30; % Nombre de simulations a faire

dt = logspace(-2,-5,nsimul);
Omega = ones(1,nsimul); %                     TODO: Choisir des valeurs de Omega pour trouver la resonance

paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = ['simulations/',paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    error = zeros(1,nsimul);
elseif strcmp(paramstr, 'Omega')
    Emax = zeros(1,nsimul);
end

g=9.81;
L=0.1;
theta0=1e-6;
w0=sqrt(g/L);

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation

    if strcmp(paramstr, 'dt')
      t = data(end,1);
      theta = data(end,2);
      theta_th = theta0*cos(w0*t);
      error(i) = abs(theta-theta_th);
    elseif strcmp(paramstr, 'Omega')
      Emax(i) = 0; %                           TODO: Calculer le maximum de l'energie
    end
end


%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    figure
    loglog(dt, error, 'k+')
    xlabel('\Delta t')
    ylabel('Erreur sur \theta(t_{fin}) [rad]')
    grid on
    print(['figures/etudeConvDt'], '-depsc');
elseif strcmp(paramstr, 'Omega')
    figure
    plot(Omega, Emax, 'k-+')
    xlabel('\Omega [rad/s]')
    ylabel('max(E_{mec}(t)) [J]')
    grid on
    print(['figures/rechercheOmega',schema], '-depsc');
end
