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

nsimul = 20; % Nombre de simulations a faire

dt = logspace(-2,-3,nsimul);
Omega = ones(1,nsimul); %                     TODO: Choisir des valeurs de Omega pour trouver la resonance
theta0 = linspace(1e-7,pi-1e-7,nsimul);

paramstr = 'theta0'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = theta0; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

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
elseif strcmp(paramstr, 'theta0')
    T = zeros(1, nsimul); 
    error = zeros(1,nsimul);
end

g=9.81;
L=0.1;
theta0petit=1e-6;
w0=sqrt(g/L);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paramstr, 'dt')
      t = data(end,1);
      theta = data(end,2);
      theta_th = theta0petit*cos(w0*t);
      error(i) = abs(theta-theta_th);
    elseif strcmp(paramstr, 'Omega')
      Emax(i) = 0; %                           TODO: Calculer le maximum de l'energie
    elseif strcmp(paramstr, 'theta0')
      t = data(:,1);
      theta = data(:,2);
      l=1;k=0;
      t_P=zeros(3);
      while (l < size(theta,1)-1 & k<3)
          if sign(theta(l))~=sign(theta(l+1))
              k=k+1;
              t_P(k)=t(l+1);
          end
          l=l+1;
      end
      T(i)=t_P(3)-t_P(1); 
      error(i)=abs(T(i)-4/w0*ellipticK(sin(theta(1)/2)^2));
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
    print(['figures/rechercheOmega'], '-depsc');
elseif strcmp(paramstr, 'theta0')
    figure
    plot(theta0, T, 'k+')
    xlabel('\theta_0 [rad]')
    ylabel('T [s]')
    grid on
    print(['figures/theta0'], '-depsc');
    
    figure
    plot(theta0, error, 'k+')
    xlabel('\theta_0 [rad]')
    ylabel('Erreur sur T [s]')
    grid on
    print(['figures/theta0error'], '-depsc');
end
