% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=0;
tFin=2*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             1         % dt
             1e-5      % precision
             "true"    % adaptatif
             "a.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);



v0=1200;
r0=314159000;
mT=5.972e24;
h=10000;
RT=6378100;
vMax_th=sqrt(v0^2+2*G*mT*(1/(h+RT)-1/r0));
alpha = pi-asin(vMax_th*(h+RT)/(v0*r0));
vx0=v0*cos(alpha);
vy0=v0*sin(alpha);

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [0               384748000       r0              % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             mT              7.3477e-22      5809            % m
             RT              1737000         1.95            % R
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);


%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable

nsimul = 10; % Nombre de simulations à faire

theta     = linspace(0,2*pi,nsimul+1);
dt        = logspace(1.5,0,nsimul); % Valeurs du parametre a scanner
precision = logspace(-1,-8,nsimul); % Valeurs du parametre a scanner

paraName='precision'; % Nom du parametre a scanner

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
    param=precision;
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

if strcmp(paraName, 'dt') || strcmp(paraName, 'precision')
    hmin = zeros(1,nsimul);
    vmax = zeros(1,nsimul);
    nsteps = ones(1,nsimul);
elseif strcmp(paramstr, 'theta')
    error = zeros(1,nsimul);
end

% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3


for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paraName, 'dt') || strcmp(paraName, 'precision')
        t = data(:,1);
        Pt = data(:,4);
        xT = data(:,5);
        yT = data(:,6);
        xA = data(:,17);
        yA = data(:,18);
        vx= data(:,20);
        vy= data(:,21);
        nsteps(i)=size(t,1)-1;
        R = T1.Terre(8);
        hmin(i)=inter_min(t,sqrt(xA.^2+yA.^2),3);
        vmax(i)=inter_max(t,sqrt(vx.^2+vy.^2),3);
%         hmin(i)=min(sqrt(xA.^2+yA.^2));
%         vmax(i)=max(sqrt(vx.^2+vy.^2));
        if strcmp(paraName, 'precision')
            dt=t(2:end-1)-t(1:end-2);
        end
    elseif strcmp(paramstr, 'theta')
        t = data(:,1);
    end
end


%% Figures %%
%%%%%%%%%%%%%
%

if strcmp(paraName, 'dt') || strcmp(paraName, 'precision')
    figure
    angle=linspace(0,2*pi);
    plot(RT*cos(angle),RT*sin(angle),'r')
    hold on
    plot(xT,yT,'+r',xA,yA)
    hold off
    axis equal
    grid on

    figure
    loglog(nsteps,abs(hmin-h-RT),'+',nsteps,1e11*nsteps.^(-4),'--')
    grid on
    xlabel('\Deltat [s]')
    ylabel('Erreur sur h_{min} [m]')

    figure
    loglog(nsteps,abs(vmax-vMax_th),'+',nsteps,2e8*nsteps.^(-4),'--')
    grid on
    xlabel('\Deltat [s]')
    ylabel('Erreur sur v_{max} [m]')

    if strcmp(paraName, 'precision')
        figure
        plot(t(1:end-2),dt)
    end
end
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

clear all;
