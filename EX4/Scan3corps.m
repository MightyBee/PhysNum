% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


rL=384748000;
G=6.674e-11;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

alpha = 0.3 %0.2138; %pi-asin(vMax_th*(h+RT)/(v0A*r0));
v0A=1200;
vx0A=v0A*cos(alpha);
vy0A=v0A*sin(alpha)+omega*rA;

vy0T=omega*rT
vy0L=omega*rL

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [rT              rL              rA              % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               vx0A            % vx0
             vy0T            vy0L            vy0A            % vy0
             0               0               0               % vz0
             mT              mL              5809            % m
             RT              RL              1.95            % R
             0               0               0            ]; % Cx

T=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T);

change_config(0,'tFin',10*24*60*60);
change_config(0,'precision',0.00001);
change_config(0,'adaptatif','true');


%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 

nsimul = 360; % Nombre de simulations à faire
paraName='theta'; % Nom du parametre a scanner

if strcmp(paraName,'theta')   
    paramstr = {"vx0"; "vy0"};
    theta = linspace(0,2*pi,nsimul+1);
    param = [v0A*cos(theta(1:nsimul)); v0A*sin(theta(1:nsimul))+omega*rA]; % Valeurs du parametre a scanner
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

if strcmp(paramstr, 'dt')
    error = zeros(1,nsimul);
elseif strcmp(paramstr, 'precision')
    Emax = zeros(1,nsimul);
elseif strcmp(paraName, 'theta')
    t  = cell(1,nsimul);
    xT = cell(1,nsimul);
    yT = cell(1,nsimul);
    xL = cell(1,nsimul);
    yL = cell(1,nsimul);
    xA = cell(1,nsimul);
    yA = cell(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paramstr, 'dt')
        t = data(end,1);
    elseif strcmp(paramstr, 'precision')
        Emec = data(:,4);
    elseif strcmp(paraName, 'theta')
        t{i}  = data(:,1);
        xT{i} = data(:,3);
        yT{i} = data(:,4);
        xL{i} = data(:,9);
        yL{i} = data(:,10);
        xA{i} = data(:,15);
        yA{i} = data(:,16);
        if(min(sqrt((xA{i}-xT{i}).^2+(yA{i}-yT{i}).^2))<=RT || min(sqrt((xA{i}-xL{i}).^2+(yA{i}-yL{i}).^2))<=RL)
            iFin=1;
            while sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>RT && sqrt((xA{i}(iFin)-xL{i}(iFin)).^2+(yA{i}(iFin)-yL{i}(iFin)).^2)>RL
                iFin=iFin+1;
            end
            t{i}  = t{i}(1:iFin);
            xT{i} = xT{i}(1:iFin);
            yT{i} = yT{i}(1:iFin);
            xL{i} = xL{i}(1:iFin);
            yL{i} = yL{i}(1:iFin);
            xA{i} = xA{i}(1:iFin);
            yA{i} = yA{i}(1:iFin);
        else
            t{i}  = 0;
            xT{i} = 0;
            yT{i} = 0;
            xL{i} = 0;
            yL{i} = 0;
            xA{i} = 0;
            yA{i} = 0;
        end
    end
end


%% Figures %%
%%%%%%%%%%%%%
% 
figure
angle=linspace(0,2*pi);
plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
hold on
for i =1:nsimul
%     plot(xT{i},yT{i},xL{i},yL{i},xA{i},yA{i})
%     hold on
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