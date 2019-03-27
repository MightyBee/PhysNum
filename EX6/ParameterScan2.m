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
methode='T'

nsimul = 20; % Nombre de simulations a faire

N = round(logspace(1,3, nsimul));

paramstr = 'N'; % Nom du parametre a scanner
param = N; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [dossier paramstr '=' num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s trivial=false methode=%s b=0.02 N1=%g N2=%g output=%s', repertoire, executable, input, methode, param(i), 5*param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

int=zeros(nsimul,1);

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load([output{i} '_phi.out']);
    phi = data(:,2);
    int(i)=phi(N(i)+1);
end

x_T=0.12./(6*N);
y_T=int;

[a,erra,yFit_T]=fit((x_T').^2,y_T);

%% G2 %%

methode='G2';

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = [dossier paramstr '=' num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s trivial=false methode=%s b=0.02 N1=%g N2=%g output=%s', repertoire, executable, input, methode, param(i), 5*param(i), output{i});
    disp(cmd)
    system(cmd);
end


for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load([output{i} '_phi.out']);
    phi = data(:,2);
    int(i)=phi(N(i)+1);
end

x_G=0.12./(6*N);
y_G=int;

[a,erra,yFit_G]=fit((x_G').^2,y_G);

%% Figures %%
%%%%%%%%%%%%%

    fig1=figure('Position',[50,50,600,450]);
    plot(x_G.^2,y_G,'b+',x_T.^2,y_T,'r+','MarkerSize',11);
    hold on;
    plot(x_G.^2,yFit_G,'b--',x_T.^2,yFit_T,'r--','MarkerSize',11);
    xlabel('h^2 [m^{2}]','FontSize', 20)
    ylabel('\phi(b) [V]','FontSize', 20)
    set(gca,'FontSize',20)
    grid on
    lgd=legend('Méthode Gauss','Méthode Trapèzes');
    set(lgd,'fontsize',14,'Location','southwest');
    print(fig1,['figures/nontrivial_conv2'], '-depsc');
