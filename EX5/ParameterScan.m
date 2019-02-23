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
executable = 'Exercice5'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

nsimul = 8; % Nombre de simulations a faire

dt = logspace(-5,-3, nsimul);

paramstr = 'dt'; % Nom du parametre a scanner
param = dt; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dt'))
    xp = 0;
    yp = 0;
    Tp = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'dt'))
        data = load([output{i} '_T.out']);
        % TODO: interpoler la temperature en (xp,yp)
        %pour trouver dans quelle cellule se trouve (xp,yp) 
        a = fix(xp/h);
        b = fix(yp/h);
        %interpolation lin�aire
        T1= data(a*N+b,3);
        T2= data((a+1)*N+b,3);
        T3= data((a+1)*N+b+1,3);
        T4= data(a*N+b+1,3);
        
        Tp(i) = ((T2-T1)-(T3-T4))/(h*h)*(xp-a*h)*(yp-b*h);
    end
end

%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paramstr,'dt'))
    figure
    plot(dt,Tp,'k+')
    xlabel('\Delta t [s]')
    ylabel(sprintf('T(%0.2f,%0.2f) [°C]',xp,yp))
    grid on
end



