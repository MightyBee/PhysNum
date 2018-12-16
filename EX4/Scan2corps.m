% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=0;
tFin=2e7;%100*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             1000        % dt
             1e-5      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);



rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

vy0T=omega*rT;
vy0L=omega*rL;

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',  'Lune',  'Apollo13' };
variables = [ rT        rL       1e10      % x0
              0         0        0         % y0
              0         0        0         % z0
              0         0        1200      % vx0
              vy0T      vy0L     0         % vy0
              0         0        0         % vz0
              mT        mL       1         % m
              RT        RL       1.95      % R
              0         0        0      ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);


%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable

nsimul = 3; % Nombre de simulations à faire
precision = logspace(-1,-5,nsimul); % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for k = 1:nsimul
    output{k} = sprintf('simulations/deuxCorps_precision=%.15g.out',precision(k));
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s configuration0.in 2 precision=%.15g output=%s', repertoire, executable, precision(k), output{k});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%

t    = cell(1,nsimul);
emec = cell(1,nsimul);
xT   = cell(1,nsimul);
yT   = cell(1,nsimul);
vxT  = cell(1,nsimul);
vyT  = cell(1,nsimul);
xL   = cell(1,nsimul);
yL   = cell(1,nsimul);
vxL  = cell(1,nsimul);
vyL  = cell(1,nsimul);
d    = cell(1,nsimul);
p    = cell(1,nsimul);
emecMAX=zeros(1,nsimul);
nsteps=zeros(1,nsimul);

% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3
for i=1:nsimul
    data=load(output{i});
    t{i}   = data(:,1);
    xT{i}  = data(:,5);
    yT{i}  = data(:,6);
    vxT{i} = data(:,8);
    vyT{i} = data(:,9);
    xL{i}  = data(:,11);
    yL{i}  = data(:,12);
    vxL{i} = data(:,14);
    vyL{i} = data(:,15);
    emec{i}= data(:,3)+0.5*mT*(vxT{i}.^2+vyT{i}.^2)+0.5*mL*(vxL{i}.^2+vyL{i}.^2);
    emecMAX(i)=max(abs(emec{i}));
    d{i}   = sqrt((xT{i}-xL{i}).^2+(yT{i}-yL{i}).^2);
    p{i}   = mT*sqrt(vxT{i}.^2+vyT{i}.^2)+mL*sqrt(vxL{i}.^2+vyL{i}.^2);
    nsteps(i)=size(t{i},1)-1;
end
%% Figures %%
%%%%%%%%%%%%%

fig1=figure('Position',[50,50,600,400]);
plot(xT{nsimul},yT{nsimul},xL{nsimul},yL{nsimul})
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
set(gca,'fontsize',15);
lgd=legend('Terre','Lune');
set(lgd,'fontsize',14,'Location','southwest');
print(fig1,'figures/deuxCorps_tajectoire', '-depsc');

fig2=figure('Position',[50,50,600,400]);
for i = 1:nsimul
    plot(t{i},d{i})
    hold on
end
hold off
grid on
xlabel('Temps [s]')
ylabel('Distance Terre-Lune [m]')
set(gca,'fontsize',15);
lgd=legend(sprintf('N_{steps}=%d',nsteps(1)),sprintf('N_{steps}=%d',nsteps(2)),sprintf('N_{steps}=%d',nsteps(3)));
set(lgd,'fontsize',14,'Location','southwest');
print(fig2,'figures/deuxCorps_d', '-depsc');

fig3=figure('Position',[50,50,600,400]);
for i = 1:nsimul
    plot(t{i},p{i})
    hold on
end
hold off
grid on
xlabel('Temps [s]')
ylabel('Quantité de mouvement [kgm/s]')
set(gca,'fontsize',15);
lgd=legend(sprintf('N_{steps}=%d',nsteps(1)),sprintf('N_{steps}=%d',nsteps(2)),sprintf('N_{steps}=%d',nsteps(3)));
set(lgd,'fontsize',14,'Location','northwest');
print(fig3,'figures/deuxCorps_p', '-depsc');

fig4=figure('Position',[50,50,600,400]);
for i = 1:nsimul
    plot(t{i},emec{i})
    hold on
end
hold off
grid on
xlabel('Temps [s]')
ylabel('Énergie mecanique [J]')
set(gca,'fontsize',15);
lgd=legend(sprintf('N_{steps}=%d',nsteps(1)),sprintf('N_{steps}=%d',nsteps(2)),sprintf('N_{steps}=%d',nsteps(3)));
set(lgd,'fontsize',14,'Location','southwest');
print(fig4,'figures/deuxCorps_emec', '-depsc');

fig5=figure('Position',[50,50,600,400]);
plot(nsteps.^(-4),emecMAX,'+',[0 12e-14],[abs(emec{nsimul}(1)) abs(emec{nsimul}(1))],'--')
grid on
xlabel('1/N_{steps}^4')
ylabel('max|E_{mec}| [J]')
set(gca,'fontsize',15);
print(fig5,'figures/deuxCorps_convEmec', '-depsc');

clear all;