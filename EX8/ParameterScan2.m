% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant deux des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur de deux parametres du fichier d'input
% par les valeurs scannees.


%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice8'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'classique.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nx = 11;
nt =11;
nsimul=nx*nt;

m=1;
omega=0.003;
P0=2*pi*14/400;
E0=P0^2/2; %mean(data(:,4));


Ninters=round(logspace(1.9,3.1,nx));
N=round(logspace(2.3,4.2, nt));
dt=5000./N;
[X,T]=meshgrid(Ninters, dt);
paramstr = {"Ninters"; "dt"}; % Nom du parametre a scanner
param = [reshape(X,[],1) reshape(T,[],1)]; % Valeurs du parametre a scanner


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(nx, nt); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nx
    for j = 1:nt
        % Sring des paramètres à varier
        parameter = '';
        for k=1:size(paramstr,1)
            parameter=[parameter sprintf('%s=%.15g ', paramstr{k}, param(i+(j-1)*nt,k))];
        end
        parameter=strip(parameter);
        % Nom du fichier de sortie
        output{i,j} = [dossier strrep(parameter, ' ', '_')];
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, input, parameter, output{i,j});
        disp(cmd)
%         system(cmd);
    end
end

%% Analyse %%
%%%%%%%%%%%%%

xmoy  = zeros(nx,nt);
pmoy  = zeros(nx,nt);
dxmoy = zeros(nx,nt);
dpmoy = zeros(nx,nt);
errx  = zeros(nx,nt);
errp  = zeros(nx,nt);

for i = 1:nx % Parcours des resultats de toutes les simulations
    for j = 1:nt
        data = load([output{i,j} '_obs.out']);
        t=data(:,1);
        x=data(:,5);
        E0=mean(data(:,4));
        xmoy(i,j) = x(end);
        x2moy   = data(end,6);
        dxmoy(i,j) = sqrt(x2moy-xmoy(i)^2);
        p=data(:,7);
        pmoy(i,j) = p(end);
        p2moy   = data(end,8);
        dpmoy(i,j) = sqrt(p2moy-pmoy(i)^2);
        A=sqrt(2*E0)/omega;
        errx(i,j)=max(abs(x/A-sin(omega*t)));
        errp(i,j)=max(abs(p/(A*omega)-cos(omega*t)));
    end
end


%% Figures %%
%%%%%%%%%%%%%

n=1;
fig(n)=figure;
surf((400./X),T,errx)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
xlabel('$\Delta h$','Interpreter','Latex','FontSize',22)
ylabel('$\Delta t$','Interpreter','Latex','FontSize',22)
zlabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0}|}$','Interpreter','Latex','FontSize',22)
xticks([0.5 1 2 3 4 5])
yticks([1 2 5 10])
zticks([0.01 0.1 1])
grid on
print(fig(n),'figures/conv_2d_x', '-depsc');
n=n+1;

fig(n)=figure;
surf((400./X),T,errp)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
set(gca, 'ZScale', 'log')
xlabel('$\Delta h$','Interpreter','Latex','FontSize',22)
ylabel('$\Delta t$','Interpreter','Latex','FontSize',22)
zlabel('$\max_t{|\frac{\langle p \rangle - p_{\rm th}}{A_0 \omega}|}$','Interpreter','Latex','FontSize',22)
xticks([0.5 1 2 3 4 5])
yticks([1 2 5 10])
zticks([0.01 0.1 1])
grid on
print(fig(n),'figures/conv_2d_p', '-depsc');
n=n+1;

% fig(n)=figure;
% hold on;
% for i=1:nt
%     plot(400./N,errx(i,:))
% end
% hold off
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% xlabel('$\Delta h$','Interpreter','Latex','FontSize',22)
% ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0 }|}$','Interpreter','Latex','FontSize',22)
% grid on, box on
% print(fig(n),'figures/conv_2d_p', '-depsc');
% n=n+1;
% 
% fig(n)=figure;
% hold on;
% for i=1:nt
%     plot(dt,errx(:,i))
% end
% hold off
% set(gca, 'XScale', 'log')
% set(gca, 'YScale', 'log')
% xlabel('$\Delta t$','Interpreter','Latex','FontSize',22)
% ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0 }|}$','Interpreter','Latex','FontSize',22)
% grid on, box on
% print(fig(n),'figures/conv_2d_p', '-depsc');
% n=n+1;