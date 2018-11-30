% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%

v=1200;
alpha=pi-10.86*pi/180;
vx0=v*cos(alpha);
vy0=v*sin(alpha);

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [0               384748000       314159000       % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             5.972e24        7.3477e8        5809            % m
             6378100         1737000         1.95            % R
             0               0               0            ]; % Cx

T=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T);

%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 

nsimul = 1; % Nombre de simulations à faire
paraName='dt'; % Nom du parametre a scanner

theta = linspace(0,2*pi,nsimul+1);
dt=logspace(8,4,nsimul); % Valeurs du parametre a scanner

if strcmp(paraName,'theta')   
    paramstr = {"vx0"; "vy0"};
    v=10;
    param = [v*cos(theta(1:nsimul)); v*sin(theta(1:nsimul))]; % Valeurs du parametre a scanner
    configfileNb=3;
elseif strcmp(paraName,'dt')
    change_config(0,'adaptatif','false');
    paramstr={"dt"};
    param=dt;
    configfileNb=0;
elseif strcmp(paraName,'precision')
    change_config(0,'adaptatif','true');
    paramstr={"precision"};
    param=logspace(-4,-8,nsimul); % Valeurs du parametre a scanner
    configfileNb=0;
end

%% Simulations %%
%%%%%%%%%%%%%%%%%




input=sprintf('configuration%d.in', configfileNb);
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for k = 1:nsimul
    parameter = "";
    for i=1:size(paramstr,1)
      parameter=parameter+sprintf('%s=%.15g ', paramstr{i,1}, param(i,k));
    end
    parameter=strip(parameter);
    output{1,k} = "simulations/"+strrep(parameter, ' ', '_')+".out";
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %d %s configuration0.in 1 output=%s', repertoire, executable, input, size(param,1), parameter, output{1,k});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

if strcmp(paraName, 'dt')
    hmin = zeros(1,nsimul);
elseif strcmp(paramstr, 'precision')
    Emax = zeros(1,nsimul);
elseif strcmp(paramstr, 'theta')
    error = zeros(1,nsimul);
end

% 1   2  3  4  5   6   7     8  9  10  11  12  13    14 15 16  17  18  19
% t   x1 y1 z1 vx1 vy1 vz1   x2 y2 z2  vx2 vy2 vz2   x3 y3 z3  vx3 vy3 vz3


for i = 1:nsimul % Parcours des resultats de toutes les simulations 
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paraName, 'dt')
        t = data(:,1);
        x = data(:,2);
        y = data(:,3);
        R = T.Terre(8);
        hmin(i)=min(sqrt(x.^2+y.^2));
    elseif strcmp(paramstr, 'precision')
        Emec = data(:,4);
    elseif strcmp(paramstr, 'theta')
        t = data(:,1);
    end
end


%% Figures %%
%%%%%%%%%%%%%
% 

figure
plot(x,y)
figure
semilogx(dt,hmin);
% if strcmp(paramstr, 'dt')
%     figure('Position',[50,50,600,400]);
% %     loglog(dt, error, 'k+')
% %     xlabel('\Deltat [s]')
% %     ylabel('Erreur sur \theta(t_{fin}) [rad]')
%     plot(dt.*dt,error,'k+')
%     xlabel('(\Deltat)^2 [s^2]')
%     ylabel('\theta(t_{fin}) [rad]')
%     set(gca,'fontsize',15);
%     title('$\Omega$=1$\omega_0$  $d$=0.04  $\kappa$=0', 'Fontweight','normal','Interpreter','latex');
%     grid on
%     print('figures/etudeConvDt', '-depsc');
% elseif strcmp(paramstr, 'Omega')
%     figure('Position',[50,50,600,400]);
%     plot(Omega, Emax, 'k+')
%     xlabel('\Omega [rad/s]')
%     ylabel('max(E_{mec}(t)) [J]')
%     set(gca,'fontsize',15);
%     grid on
%     print('figures/rechercheOmega', '-depsc');
% elseif strcmp(paramstr, 'theta0')
%     fig1=figure('Position',[50,50,600,400]);
%     plot(theta0_ana, T_ana,'r-',theta0, T_num, 'k+')
%     lgd=legend('Analytique','Numérique');
%     set(lgd,'fontsize',14,'Location','northwest');
%     xlabel('\theta_0 [rad]')
%     ylabel('T [s]')
%     set(gca,'fontsize',15);
%     grid on
%     print(fig1,'figures/theta0', '-depsc');
% 
%     fig2=figure('Position',[50,50,600,400]);
%     plot(theta0, error, 'k+')
%     xlabel('\theta_0 [rad]')
%     ylabel('Erreur sur T [s]')
%     set(gca,'fontsize',15);
%     grid on
%     print(fig2,'figures/theta0error', '-depsc');
% end
