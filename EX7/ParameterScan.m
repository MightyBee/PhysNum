% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.


%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 20; % Nombre de simulations a faire


L=20;
u=6;
% N = round((logspace(1,3, nsimul)./4))*4+1;
N=round(logspace(2,3, nsimul));


Omega=linspace(0,5,nsimul);
Tfin=linspace(20,150,nsimul);
paramstr = 'tfin'; % Nom du parametre a scanner
param = Tfin; % Valeurs du parametre a scanner


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
        f = data(:,2:end);
%         y(i) = f(end,(N(i)-1)/4+1);
        a = floor(x0/dx(i))
        y(i)=f(end,a)+(f(end,a+1)-f(end,a))/dx(i)*(x0-dx(i)*a);
    end
end

% [a,erra,yFit]=fit(dt',Tp');
if(strcmp(paramstr,'omega'))
      y = zeros(1,nsimul);
    for i=1:nsimul
        data =  load([output{i} '_E.out']);  
       y(i)=max(data(:,2));
      % y(i)=max(E);
      %  for j=1:7961
          %  if data(j,2)>y(i)
            %    y(i)=data(j,2)
           % end
       % end
    end
end

if(strcmp(paramstr,'tfin'))
      y = zeros(1,nsimul);
    for i=1:nsimul
        data =  load([output{i} '_E.out']);  
       y(i)=max(data(:,2));
      % y(i)=max(E);
      %  for j=1:7961
          %  if data(j,2)>y(i)
            %    y(i)=data(j,2)
           % end
       % end
    end
end
%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paramstr,'Npoints'))
    figure
    h=plot(1./N,y,'k+');
    xlabel('$\Delta x \ \rm [m]$','Interpreter','Latex','FontSize', 20)
    ylabel(sprintf('T[]',xp,yp),'FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
    grid on
    lgd=legend('Valeurs numériques', 'Régression linéaire');
    set(lgd,'fontsize',14,'Location','northwest');
%     print(fig1,'conv', '-depsc');
end

if(strcmp(paramstr,'omega'))
    figure
    h=plot(Omega,y,'k+');
    grid on
end

if(strcmp(paramstr,'tfin'))
    figure
    h=plot(Tfin,y,'k+');
    grid on
end