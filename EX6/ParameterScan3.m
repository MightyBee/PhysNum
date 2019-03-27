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
    cmd = sprintf('%s%s %s methode=G3 trivial=false b=0.02 N1=%g N2=%g output=%s', repertoire, executable, input, param(i), 5*param(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

int=zeros(nsimul,1);
tri=zeros(nsimul,1);
rho=cell(nsimul,1);
rcaca=cell(nsimul,1);
for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load([output{i} '_Er_Dr.out']);
    rmid = data(:,1);
    Er = data(:,2);
    Dr = data(:,3);
    data = load([output{i} '_phi.out']);
    r = data(:,1);
    phi = data(:,2);
    data = load([output{i} '_rholib_divEr_divDr.out']);
    rmidmid = data(:,1);
    rholib = data(:,2);
    divEr = data(:,3);
    divDr = data(:,4);
    rho{i}=divEr-divDr;
    rcaca{i}=rmidmid;
    int(i)=rho{i}(N(i));
    tri(i)=(rho{i}(N(i))+rho{i}(N(i)-1))*(rmidmid(N(i))-rmidmid(N(i)-1))/2+(rho{i}(N(i)+1)+rho{i}(N(i)))*(rmidmid(N(i)+1)-rmidmid(N(i)))/2;
    fprintf('%.15g \n',rmidmid(N(i)));
end

x=0.02./N;
y=int;


[a,erra,yFit]=fit(x',tri);

%% Figures %%
%%%%%%%%%%%%%
fig4=figure('Position',[50,50,600,450]);
hold on
for i=1:nsimul
plot(rcaca{i},rho{i},'DisplayName',sprintf('N_1 + N_2 = %d',6*N(i)))
end
xlabel('r')
ylabel('\rho/\epsilon_0 [V/m^2]')
grid on, box on
set(gca,'FontSize',20)
lgd=legend('show');
set(lgd,'fontsize',14,'Location','southeast');
%print(fig4,['figures/rhcacao_' filename], '-depsc');

    fig1=figure('Position',[50,50,600,450]);
    h=plot(x,y,'+');
    xlabel('h [m]','FontSize', 20)
    ylabel('\rho_{pol}(b) [V]','FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
    grid on
    %lgd=legend('Valeurs numériques', 'Régression linéaire');
    %set(lgd,'fontsize',14,'Location','southwest');
    print(fig1,'figures/nontrivial_divergence_rhoLib', '-depsc');
fig1=figure('Position',[50,50,600,450]);
    h=plot(x,tri,'+')%,x,yFit,'--')%,x.^2,yFit,'--');
    xlabel('$h \ \rm [m]$','Interpreter','Latex','FontSize', 20)
    ylabel('$\sigma_{pol}(b) \ \rm [V]$','Interpreter','Latex','FontSize', 20)
    set(gca,'FontSize',20)
    set(h,'MarkerSize',11)
    grid on
    %lgd=legend('Valeurs numériques', 'Régression linéaire');
    %set(lgd,'fontsize',14,'Location','southwest');
    print(fig1,'figures/nontrivial_charge_b', '-depsc');
