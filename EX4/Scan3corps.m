% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


% Paramètres généraux %%%%

G=6.674e-11;
rho0=1.2;
tFin=20*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             50        % dt
             1e-5      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);


% Paramètres des corps %

rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

v0A=1200;
vy0T=omega*rT
vy0L=omega*rL

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [rT              rL              rA              % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               v0A             % vx0
             vy0T            vy0L            0               % vy0
             0               0               0               % vz0
             mT              mL              5809            % m
             RT              RL              1.95            % R
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);



%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 


nsimul = 360; % Nombre de simulations à faire


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraName='theta'; % Nom du parametre a scanner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trajectoire=true;


if strcmp(paraName,'theta')   
    paramstr = {"vx0"; "vy0"};
    theta = linspace(0,2*pi,nsimul+1); %[linspace(3.676,3.680,round(nsimul/3)) linspace(4.263,4.267,round(nsimul/3)) linspace(5.062,5.066,nsimul-2*round(nsimul/3))];
    theta = theta(1:nsimul);
    param = [v0A*cos(theta); v0A*sin(theta)+omega*rA]; % Valeurs du parametre a scanner
    configfileNb=3;
elseif strcmp(paraName,'condIn')   
    paramstr = {"vx0"; "vy0"};
    nsimul=round(sqrt(nsimul));
    theta = linspace(0,2*pi,nsimul+1);
    theta = theta(1:nsimul);
    v0A = linspace(10,10000,nsimul);
    param=[];
    for i=1:nsimul
       param = [[param] [v0A(i)*cos(theta); v0A(i)*sin(theta)+omega*rA]]; % Valeurs du parametre a scanner
    end
    nsimul = nsimul^2;
    configfileNb=3;
elseif strcmp(paraName,'dt')
    change_config(0,'adaptatif','false');
    paramstr={"dt"};
    param=logspace(8,4,nsimul); % Valeurs du parametre a scanner
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
    % Sring des paramètres à varier
    parameter = "";
    for i=1:size(paramstr,1)
      parameter=parameter+sprintf('%s=%.15g ', paramstr{i,1}, param(i,k));
    end
    parameter=strip(parameter);
    % Nom du fichier de sortie
    output{1,k} = "simulations/"+strrep(parameter, ' ', '_')+".out";
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %d %s configuration0.in 1 output=%s', repertoire, executable, input, size(param,1), parameter, output{1,k});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

if strcmp(paramstr, 'dt')
    error = zeros(1,nsimul);
elseif strcmp(paramstr, 'precision')
    Emax = zeros(1,nsimul);
elseif strcmp(paraName, 'theta')
    t  = cell(1,nsimul);
    acc= zeros(1,nsimul);
    Pt = cell(1,nsimul);
    xT = cell(1,nsimul);
    yT = cell(1,nsimul);
    xL = cell(1,nsimul);
    yL = cell(1,nsimul);
    xA = cell(1,nsimul);
    yA = cell(1,nsimul);
    hmin= zeros(1,nsimul);
    lmin= zeros(1,nsimul);
elseif strcmp(paraName, 'condIn')
    t  = cell(1,nsimul);
    acc= zeros(sqrt(nsimul));
    Pt = cell(1,nsimul);
    xT = cell(1,nsimul);
    yT = cell(1,nsimul);
    xL = cell(1,nsimul);
    yL = cell(1,nsimul);
    xA = cell(1,nsimul);
    yA = cell(1,nsimul);
end



% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paramstr, 'dt')
        t = data(end,1);
    elseif strcmp(paramstr, 'precision')
        Emec = data(:,3);
    elseif strcmp(paraName, 'theta') 
        t{i}  = data(:,1);
%         acc(i)= inter_max(t{i},data(:,2),3);
        acc(i)= max(data(:,2));
        Pt{i} = data(:,4);
        xT{i} = data(:,5);
        yT{i} = data(:,6);
        xL{i} = data(:,11);
        yL{i} = data(:,12);
        xA{i} = data(:,17);
        yA{i} = data(:,18);
        hmin(i)=min(sqrt((xA{i}-xT{i}).^2+(yA{i}-yT{i}).^2));
        lmin(i)=min(sqrt((xA{i}-xL{i}).^2+(yA{i}-yL{i}).^2));
        if hmin(i)<=RT && lmin(i)>RL
            iFin=1;
            while sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>RT %&& sqrt((xA{i}(iFin)-xL{i}(iFin)).^2+(yA{i}(iFin)-yL{i}(iFin)).^2)>RL
                iFin=iFin+1;
            end
            t{i}  = t{i}(1:iFin);
%             acc(i)= inter_max(t{i},data(1:iFin,2),3);
            acc(i)= max(data(1:iFin,2));
            Pt{i} = Pt{i}(1:iFin);
            xT{i} = xT{i}(1:iFin);
            yT{i} = yT{i}(1:iFin);
            xL{i} = xL{i}(1:iFin);
            yL{i} = yL{i}(1:iFin);
            xA{i} = xA{i}(1:iFin);
            yA{i} = yA{i}(1:iFin);
        else
            t{i}  = nan;
            acc(i)= nan;
            Pt{i} = nan;
            xT{i} = nan;
            yT{i} = nan;
            xL{i} = nan;
            yL{i} = nan;
            xA{i} = nan;
            yA{i} = nan;
        end
    elseif strcmp(paraName, 'condIn')
        t{i}  = data(:,1);
        acc(floor((i-1)/sqrt(nsimul))+1,mod(i-1,sqrt(nsimul))+1)= max(data(:,2));
        Pt{i} = data(:,4);
        xT{i} = data(:,5);
        yT{i} = data(:,6);
        xL{i} = data(:,11);
        yL{i} = data(:,12);
        xA{i} = data(:,17);
        yA{i} = data(:,18);
%         if min(sqrt((xA{i}-xT{i}).^2+(yA{i}-yT{i}).^2))<=RT && min(sqrt((xA{i}-xL{i}).^2+(yA{i}-yL{i}).^2))>RL
%             iFin=1;
%             while sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>RT %&& sqrt((xA{i}(iFin)-xL{i}(iFin)).^2+(yA{i}(iFin)-yL{i}(iFin)).^2)>RL
%                 iFin=iFin+1;
%             end
%             t{i}  = t{i}(1:iFin);
%             acc(floor((i-1)/sqrt(nsimul))+1,mod(i-1,sqrt(nsimul))+1)= max(data(1:iFin,2));
%             Pt{i} = Pt{i}(1:iFin);
%             xT{i} = xT{i}(1:iFin);
%             yT{i} = yT{i}(1:iFin);
%             xL{i} = xL{i}(1:iFin);
%             yL{i} = yL{i}(1:iFin);
%             xA{i} = xA{i}(1:iFin);
%             yA{i} = yA{i}(1:iFin);
%         else
%             t{i}  = 0;
%             acc(floor((i-1)/sqrt(nsimul))+1,mod(i-1,sqrt(nsimul))+1)= 0;
%             Pt{i} = 0;
%             xT{i} = 0;
%             yT{i} = 0;
%             xL{i} = 0;
%             yL{i} = 0;
%             xA{i} = 0;
%             yA{i} = 0;
%         end
    end
end


%% Figures %%
%%%%%%%%%%%%%
% 
if strcmp(paraName, 'theta')
    if trajectoire
        figure
        angle=linspace(0,2*pi,10000);
        plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
        hold on
        for i =1:nsimul
%             plot(xT{i},yT{i},xL{i},yL{i},xA{i},yA{i})
%             hold on
            plot(xT{i}.*cos(omega*t{i})+yT{i}.*sin(omega*t{i}), -xT{i}.*sin(omega*t{i})+yT{i}.*cos(omega*t{i}), 'k+')
            hold on
            plot(xL{i}.*cos(omega*t{i})+yL{i}.*sin(omega*t{i}), -xL{i}.*sin(omega*t{i})+yL{i}.*cos(omega*t{i}), 'r+')
            hold on
            plot(xA{i}.*cos(omega*t{i})+yA{i}.*sin(omega*t{i}), -xA{i}.*sin(omega*t{i})+yA{i}.*cos(omega*t{i}))
            hold on
        end
        hold off
        xlabel('x [m]')
        ylabel('y [m]')
        grid on
        axis equal
    end 

    figure
    plot(theta,acc,'p')
    xlabel('\theta_0 [rad]')
    ylabel('max(acc) [m/s^2]')
    grid on
    
    figure
    plot(theta,hmin-RT,'p')
    xlabel('\theta_0 [rad]')
    ylabel('min(h) [m/s^2]')
    grid on

    figure
    plot(t{50},Pt{50},'p')
    xlabel('Temps [s]')
    ylabel('Puissance de traînée [J]')
    grid on

    figure
    plot(t{50},sqrt(xA{50}.^2+yA{50}.^2),'p')
    xlabel('Temps [s]')
    ylabel('Distance au centre de masse [m]')
    grid on
elseif strcmp(paraName, 'condIn')
    if trajectoire
        figure
        angle=linspace(0,2*pi);
        plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
        hold on
        for i =1:nsimul
%             plot(xT{i},yT{i},xL{i},yL{i},xA{i},yA{i})
%             hold on
            plot(xT{i}.*cos(omega*t{i})+yT{i}.*sin(omega*t{i}), -xT{i}.*sin(omega*t{i})+yT{i}.*cos(omega*t{i}), 'k+')
            hold on
            plot(xL{i}.*cos(omega*t{i})+yL{i}.*sin(omega*t{i}), -xL{i}.*sin(omega*t{i})+yL{i}.*cos(omega*t{i}), 'r+')
            hold on
            plot(xA{i}.*cos(omega*t{i})+yA{i}.*sin(omega*t{i}), -xA{i}.*sin(omega*t{i})+yA{i}.*cos(omega*t{i}))
            hold on
        end
        hold off
        xlabel('x [m]')
        ylabel('y [m]')
        grid on
        axis equal
    end
    figure
    surf(theta,v0A,acc);
    xlabel('\theta_0 [rad]')
    ylabel('v_0 [m/s]')
    zlabel('max(acc) [m/s^2]')
    grid on
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