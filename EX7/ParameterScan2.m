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
input = 'configuration.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 20; % Nombre de simulations a faire


% N = round((logspace(1,3, nsimul)./4))*4+1;
cfl=linspace(0.99,1.01,nsimul);

paramstr = 'CFL'; % Nom du parametre a scanner
param = cfl; % Valeurs du parametre a scanner


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

    E = zeros(1,nsimul);


for i = 1:nsimul % Parcours des resultats de toutes les simulations
    
data = load([output{i},'_E.out']);
E(i) = max(data(:,2));
end

% [a,erra,yFit]=fit(dt',Tp');

%% Figures %%
%%%%%%%%%%%%%

    fig2=figure('Position',[50,50,600,450]);
semilogy(cfl,E,'+')
xlabel('$\beta_{\rm CFL}$','Interpreter','Latex')
ylabel('$E_{\max} \ \rm [J]$','Interpreter','Latex')
grid on
set(gca,'FontSize',25)
% pos=get(gca,'position');  % retrieve the current values
% pos(3)=0.9*pos(3);        % try reducing width 10%
% set(gca,'position',pos);  % write the new values
print(fig2,sprintf('figures/energie_scan',fichier), '-depsc');
% print(fig1,sprintf('figures/%s',fichier), '-depsc');

