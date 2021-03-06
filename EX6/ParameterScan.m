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
executable = 'Exercice6'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 12; % Nombre de simulations a faire

N = round(logspace(1,4, nsimul));

paramstr = 'N'; % Nom du parametre a scanner
param = N; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [dossier paramstr '=' num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s trivial=true b=0.06 N1=%g N2=%g output=%s', repertoire, executable, input, param(i), param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

int=zeros(nsimul,1);

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load([output{i} '_Er_Dr.out']);
    rmid = data(:,1);
    Er = data(:,2);
    Dr = data(:,3);
    data = load([output{i} '_phi.out']);
    r = data(:,1);
    phi = data(:,2);
    int(i)=phi(1);
    data = load([output{i} '_rholib_divEr_divDr.out']);
    rmidmid = data(:,1);
    rholib = data(:,2);
    divEr = data(:,3);
    divDr = data(:,4);
end

x=0.12./(2*N);
y=abs(int-0.12^2/4);

[a,erra,yFit]=fit(log(x'),log(y));
fprintf('Ordre de convergence : %.5g \n',a)

%% Figures %%
%%%%%%%%%%%%%

    fig1=figure('Position',[50,50,600,450]);
    h=loglog(x,y,'+',x,exp(yFit),'--',x,exp(2*log(x)-2),'--');
    xlabel('h [m]','FontSize', 20)
    ylabel('|\phi(0)-\phi_0| [V]','FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
    grid on
    lgd=legend('Valeurs numériques', 'Régression linéaire', '\propto h^2');
    set(lgd,'fontsize',15,'Location','southeast');
    print(fig1,'figures/trivial_conv', '-depsc');
