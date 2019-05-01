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
input = 'scan_const.in'; % Nom du fichier d'entree de base
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
    dt = zeros(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paramstr,'Npoints'))
        data = load([output{i} '_f.out']);
        t = data(end-1:end,1);
        dt(i)=abs(t(2)-t(1));
        f = data(end-1:end,2:end);
%         y(i) = f(end,(N(i)-1)/4+1);
        a = floor(x0/dx(i)+1);
        y1=f(1,a)+(f(1,a+1)-f(1,a))/dx(i)*(x0-dx(i)*a);
        y2=f(2,a)+(f(2,a+1)-f(2,a))/dx(i)*(x0-dx(i)*a);
        y(i)=y1+(y2-y1)/(t(2)-t(1))*(1.5-t(1));
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
    fig1=figure('Position',[50,50,600,450]);
    h=loglog(20./N,err,'k+');
    xlabel('$\Delta x \ \rm [m]$','Interpreter','Latex')
    ylabel('Erreur [m]','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig1,'figures/conv0_dx', '-depsc');

    fig2=figure('Position',[50,50,600,450]);
    h=loglog(dt,err,'k+');
    xlabel('$\Delta t \ \rm [s]$','Interpreter','Latex')
    ylabel('Erreur [m]','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig2,'figures/conv0_dt', '-depsc');
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
