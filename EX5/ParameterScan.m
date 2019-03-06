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
dossier='simulations/';

nsimul = 100; % Nombre de simulations a faire

<<<<<<< HEAD
dt = linspace(0.000775,0.0008, nsimul);
=======
dt = logspace(-5,-3, nsimul);
>>>>>>> 26360ee61befe250ed5205a1acdec20c827bf400

paramstr = 'dt'; % Nom du parametre a scanner
param = dt; % Valeurs du parametre a scanner

N = 80;
L=0.1;
h=L/N;
%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [dossier,paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paramstr,'dt'))
    xp = 0.06;
    yp = 0.03;
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
<<<<<<< HEAD
        T1= data(a*(N+1)+b,3);
        T2= data((a+1)*(N+1)+b,3);
        T3= data((a+1)*(N+1)+b+1,3);
        T4= data(a*(N+1)+b+1,3);
        
=======
        T1= data(a*N+b,3);
        T2= data((a+1)*N+b,3);
        T3= data((a+1)*N+b+1,3);
        T4= data(a*N+b+1,3);

>>>>>>> 26360ee61befe250ed5205a1acdec20c827bf400
        Tb=T1+(T4-T1)*(xp-a*h)/h;
        Th=T2+(T3-T2)*(xp-a*h)/h;
        Tp(i) = Tb + (Th-Tb)*(yp-b*h)/h ;
    end
end

[a,erra,yFit]=fit(dt',Tp');


%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paramstr,'dt'))
<<<<<<< HEAD
    figure
    h=plot(dt,Tp,'k+')
    xlabel('\Delta t [s]','FontSize', 20)
    ylabel(sprintf('T[�C]',xp,yp),'FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
=======
    fig1=figure('Position',[50,50,600,400]);
    plot(dt,Tp,'k+',dt,yFit,'--')
    xlabel('\Delta t [s]')
    ylabel(sprintf('T(%0.2f,%0.2f) [°C]',xp,yp))
    set(gca,'fontsize',13);
>>>>>>> 26360ee61befe250ed5205a1acdec20c827bf400
    grid on
    lgd=legend('Valeurs numériques', 'Régression linéaire');
    set(lgd,'fontsize',14,'Location','northwest');
    print(fig1,'conv', '-depsc');
end
