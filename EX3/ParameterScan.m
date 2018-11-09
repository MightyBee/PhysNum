




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

dt = logspace(-2,-5,nsimul);
Omega = linspace(9.6,10,nsimul); %                     TODO: Choisir des valeurs de Omega pour trouver la resonance
theta0 = linspace(1e-7,pi-1e-7,nsimul);

paramstr = 'Omega'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = Omega; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

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
    T_num = zeros(1, nsimul);
    theta0_ana=linspace(1e-7,pi-1e-1,200);
    T_ana=4/w0*ellipticK(sin(theta0_ana/2).*sin(theta0_ana/2));
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
<<<<<<< HEAD
        v=-10;
        for l=(data(:,4))
            if l>v

                v=l;
=======
        Emec=data(:,4);
        v=Emec(1);
        for l=1:size(Emec,1)
            if Emec(l)>v
                v=Emec(l);
>>>>>>> 30c57a5bf2cde01101cd7afee9ff0263328e1f7a
            end
        end
        Emax(i)= v; % TODO: Calculer le maximum de l'energie
    elseif strcmp(paramstr, 'theta0')
      t = data(:,1);
      theta = data(:,2);
      l=1;k=0;
      t_P=zeros(3);
      while (l < size(theta,1)-1 && k<3)
          if sign(theta(l))~=sign(theta(l+1))
              k=k+1;
              t_P(k)=t(l)+ (t(l+1)-t(l))*abs(theta(l))/abs(theta(l+1)-theta(l));
          end
          l=l+1;
      end
      T_num(i)=t_P(3)-t_P(1);
      error(i)=abs(T(1,i)-T(2,i));
    end
end


%% Figures %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    figure('Position',[50,50,600,400]);
    loglog(dt, error, 'k+')
    xlabel('\Delta t')
    ylabel('Erreur sur \theta(t_{fin}) [rad]')
    set(gca,'fontsize',15);
    grid on
    print('figures/etudeConvDt', '-depsc');
elseif strcmp(paramstr, 'Omega')
    figure('Position',[50,50,600,400]);
    plot(Omega, Emax, 'k+')
    xlabel('\Omega [rad/s]')
    ylabel('max(E_{mec}(t)) [J]')
    set(gca,'fontsize',15);
    grid on
    print('figures/rechercheOmega', '-depsc');
elseif strcmp(paramstr, 'theta0')
    fig1=figure('Position',[50,50,600,400]);
    plot(theta0_ana, T_ana,'r-',theta0, T_num, 'k+')
    lgd=legend('Analytique','Numérique');
    set(lgd,'fontsize',14,'Location','northwest');
    xlabel('\theta_0 [rad]')
    ylabel('T [s]')
    set(gca,'fontsize',15);
    grid on
    print(fig1,'figures/theta0', '-depsc');

    fig2=figure('Position',[50,50,600,400]);
    plot(theta0, error, 'k+')
    xlabel('\theta_0 [rad]')
    ylabel('Erreur sur T [s]')
    set(gca,'fontsize',15);
    grid on
    print(fig2,'figures/theta0error', '-depsc');
end
