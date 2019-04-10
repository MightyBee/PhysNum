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
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'scan_const.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 10; % Nombre de simulations a faire


% N = round((logspace(1,3, nsimul)./4))*4+1;
N=round(logspace(2,4.5, nsimul));

paramstr = 'Npoints'; % Nom du parametre a scanner
param = N; % Valeurs du parametre a scanner

L=20;
dx=L./(N-1);
x0=5;

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

if(strcmp(paramstr,'Npoints'))
    xp = 0.06;
    yp = 0.03;
    y = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'Npoints'))
        data = load([output{i} '_f.out']);
        t = data(end-1:end,1);
        f = data(end-1:end,2:end);
%         y(i) = f(end,(N(i)-1)/4+1);
        a = floor(x0/dx(i)+1);
        y1=f(1,a)+(f(1,a+1)-f(1,a))/dx(i)*(x0-dx(i)*a);
        y2=f(2,a)+(f(2,a+1)-f(2,a))/dx(i)*(x0-dx(i)*a);
        y(i)=y1+(y2-y1)/(t(2)-t(1))*(1.5-t(1));
    end
end

% [a,erra,yFit]=fit(dt',Tp');

err=abs(-sin(5/6*5-5*1.5)-y);
%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paramstr,'Npoints'))
    figure
    h=loglog(1./N,err,'k+');
    xlabel('$\Delta x \ \rm [m]$','Interpreter','Latex','FontSize', 20)
    ylabel(sprintf('T[]',xp,yp),'FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
    grid on
    lgd=legend('Valeurs numériques', 'Régression linéaire');
    set(lgd,'fontsize',14,'Location','northwest');
%     print(fig1,'conv', '-depsc');
end
