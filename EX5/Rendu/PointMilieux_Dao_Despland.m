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
dossier="simulations/";

nsimul = 20; % Nombre de simulations a faire

N = 80;
L=0.1;
h=L/N;
kappa=1.2;

xa = 2*h;
xb = xa+0.02;
xd = L-2*h;
xc = xd-0.01;


ind_i=floor((xb+xc)/2/h);
x_m=(ind_i+0.5)*h;
ind_j=floor(N/2);
y_m=(ind_j+0.5)*h;
h1=(x_m-h-xb)/(nsimul-1);
h2=(xc-x_m-h)/(nsimul-1);

d=(h1+h2)*linspace(nsimul-1,0,nsimul)+2*h;


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    parameter=sprintf('xa=%.15g xb=%.15g xc=%.15g xd=%.15g', xa+(i-1)*h1, xb+(i-1)*h1, xc-(i-1)*h2, xd-(i-1)*h2);
    output{i} = dossier+strrep(parameter, ' ', '_');
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s tfin=10 N=%d %s output=%s', repertoire, executable, input, N, parameter, output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

jx = zeros(1,nsimul);
jy = zeros(1,nsimul);


for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}+"_T.out");
    T1= data(ind_i*(N+1)+ind_j,3);
    T2= data((ind_i+1)*(N+1)+ind_j,3);
    T3= data((ind_i+1)*(N+1)+ind_j+1,3);
    T4= data(ind_i*(N+1)+ind_j+1,3);
    jx(i)=-kappa*0.5*(T2-T1+T3-T4)/h;
    jy(i)=-kappa*0.5*(T4-T1+T3-T2)/h;
end

disp(data(ind_i*(N+1)+ind_j,1:2));
disp(data((ind_i+1)*(N+1)+ind_j,1:2));
disp(data((ind_i+1)*(N+1)+ind_j+1,1:2));
disp(data(ind_i*(N+1)+ind_j+1,1:2));

[a,erra,yFit]=fit(1./d(1:round(0.8*nsimul))',jx(1:round(0.8*nsimul))');

%% Figures %%
%%%%%%%%%%%%%

jx=sqrt(jx.^2+jy.^2);

fig1=figure('Position',[50,50,500,400]);
plot(d,jx,'k+')
xlabel('d [m]')
ylabel(sprintf('|j|(%0.2f,%0.2f) [W/m]',x_m,y_m))
xlim([0 0.07])
set(gca,'fontsize',14);
grid on
print(fig1,'milieu', '-depsc');

fig2=figure('Position',[50,50,500,400]);
plot(1./d(1:round(0.8*nsimul)),jx(1:round(0.8*nsimul)),'k+',1./d(1:round(0.8*nsimul)),yFit,'--')
xlabel('d^{-1} [m^{-1}]')
ylabel(sprintf('|j|(%0.2f,%0.2f) [W/m]',x_m,y_m))
set(gca,'fontsize',14);
grid on
lgd=legend('Valeurs numériques', 'Régression linéaire');
set(lgd,'fontsize',14,'Location','northwest');
print(fig2,'milieu_inverse', '-depsc');
