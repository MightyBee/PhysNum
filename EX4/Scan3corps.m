% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


% Paramètres généraux %%%%

G=6.674e-11;
rho0=1.2;
tFin=17*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             50        % dt
             1e-6      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             1      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);


% Paramètres des corps %

rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2)

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

fprintf('%13s : DONE \n','Configuration')

%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 


nsimul = 210; % Nombre de simulations à faire


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraName='theta'; % Nom du parametre a scanner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trajectoire=true;
dim3=false;
interAcc=true;


if strcmp(paraName,'theta')   
    paramstr = {"vx0"; "vy0"};
    theta = [ linspace(0.21,0.29,nsimul/3) ; linspace(3.6,4.3,nsimul/3) ; linspace(5.0,5.7,nsimul/3)];
    theta = theta(1:nsimul);
%     theta = linspace(3.6795,3.6805,nsimul);
%     theta = linspace(4.261,4.2635,nsimul); %200 simulations (1e-6) -> indice 76, donne une double entrée atmoshere
%     theta = linspace(4.2615,4.2625,nsimul);
%    theta = linspace(5.067,5.068,nsimul);
%     theta = linspace(5.6195,5.6205,nsimul);
%     theta1= [0.18 0.32; 3.6 4.4; 5 5.7];
%     theta1 = [3.67 3.72; 4.22 4.27; 5.06 5.11; 5.58 5.63];
%     theta1= [3.86 3.68015; 4.2617 4.2618; 5.0673 5.0675; 5.62006 5.620014];
%     theta = create_linspace(theta1,nsimul);
    param = [v0A*cos(theta); v0A*sin(theta)+omega*rA]; % Valeurs du parametre a scanner
    configfileNb=3;
elseif strcmp(paraName,'condIn')   
    paramstr = {"vx0"; "vy0"};
    nsimul=round(sqrt(nsimul));
    theta = linspace(0,1/2*pi,nsimul+1);
    theta = theta(1:nsimul);
    v0A = linspace(10,10000,nsimul);
    param=[];
    for i=1:nsimul
       param = [[param] [v0A(i)*cos(theta); v0A(i)*sin(theta)+omega*rA]]; % Valeurs du parametre a scanner
    end
    nsimul = nsimul^2;
    configfileNb=3;
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

fprintf('%13s : DONE \n','Simulations')

%% Analyse %%
%%%%%%%%%%%%%

if strcmp(paraName, 'theta')
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
    vAbs= cell(1,nsimul);
    vRel= cell(1,nsimul);
    rentree=zeros(1,nsimul);
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
    if strcmp(paraName, 'theta') 
        t{i}  = data(:,1);
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
        vAbs{i}= sqrt(data(:,20).^2+data(:,21).^2);
        vRel{i}= sqrt((data(:,20)-data(:,8)).^2+(data(:,21)-data(:,9)).^2);
        if hmin(i)<=RT && lmin(i)>RL
            distAtmo=RT+120000;
            rentree(i)=1;
            iFin=1;
            while sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>distAtmo 
                iFin=iFin+1;
            end
            iAtmo=iFin;
            while sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>RT
                if sqrt((xA{i}(iFin)-xT{i}(iFin)).^2+(yA{i}(iFin)-yT{i}(iFin)).^2)>distAtmo
                    rentree(i)=0;
                    break;
                end
                iFin=iFin+1;
            end
            t{i}  = t{i}(1:iFin);
            acc(i)= inter_max(t{i},data(1:iFin,2),3);
            Pt{i} = Pt{i}(1:iFin);
            xT{i} = xT{i}(1:iFin);
            yT{i} = yT{i}(1:iFin);
            xL{i} = xL{i}(1:iFin);
            yL{i} = yL{i}(1:iFin);
            xA{i} = xA{i}(1:iFin);
            yA{i} = yA{i}(1:iFin);
            vAbs{i}= vAbs{i}(1:iFin);
            vRel{i}= vRel{i}(1:iFin);
        end
        if not (hmin(i)<=RT && lmin(i)>RL && rentree(i)==1)
            t{i}  = nan;
            acc(i)= nan;
            Pt{i} = nan;
            xT{i} = nan;
            yT{i} = nan;
            xL{i} = nan;
            yL{i} = nan;
            xA{i} = nan;
            yA{i} = nan;
            vAbs{i}= nan;
            vRel{i}= nan;
        else
            hmin(i)=nan;
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
if interAcc
    format long
    [minAcc,itt]=inter_min(theta,acc,3)
end

fprintf('%13s : DONE \n','Analyse')

%% Figures %%
%%%%%%%%%%%%%
% 
if strcmp(paraName, 'theta')
    if trajectoire
        if dim3
            figure
            for i =1:nsimul
                if t{i}(end) > 3.2e0
                    plot3(xT{i},yT{i},t{i},'r+',xL{i},yL{i},t{i},'k+',xA{i},yA{i},t{i})
                    hold on
                end
            end
            hold off
            xlabel('x [m]')
            ylabel('y [m]')
            zlabel('t [s]')
            grid on
            h = get(gca,'DataAspectRatio')
            if h(3)==1
                set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
            else
                set(gca,'DataAspectRatio',[1 1 h(3)])
            end
            
            figure
            hold on
            for i =1:nsimul
                if t{i}(end) > 3.2e0
                    plot3(xT{i}.*cos(omega*t{i})+yT{i}.*sin(omega*t{i}), -xT{i}.*sin(omega*t{i})+yT{i}.*cos(omega*t{i}),t{i},'k+')%(end)*ones(size(t{i})), 'k+')
                    hold on
                    plot3(xL{i}.*cos(omega*t{i})+yL{i}.*sin(omega*t{i}), -xL{i}.*sin(omega*t{i})+yL{i}.*cos(omega*t{i}),t{i},'r+')%(end)*ones(size(t{i})), 'r+')
                    hold on
                    plot3(xA{i}.*cos(omega*t{i})+yA{i}.*sin(omega*t{i}), -xA{i}.*sin(omega*t{i})+yA{i}.*cos(omega*t{i}),t{i})%(end)*ones(size(t{i})))
                    hold on
                end
            end
            hold off
            xlabel("x' [m]")
            ylabel("y' [m]")
            zlabel('t [s]')
            grid on
            h = get(gca,'DataAspectRatio')
            if h(3)==1
                set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))])
            else
                set(gca,'DataAspectRatio',[1 1 h(3)])
            end
        else
            fig1=figure('Position',[50,50,600,400]);
            for i =1:nsimul
                if t{i}(end) > 3.2e0
                    plot(xT{i},yT{i},'r',xL{i},yL{i},'k',xA{i},yA{i})
                    hold on
                end
            end
            hold off
            xlabel('x [m]')
            ylabel('y [m]')
            axis equal
            grid on
            set(gca,'fontsize',15);
            print(fig1,'figures/troisCorpsAdaptRho0_trajectoireR', '-depsc');
            
            fig2=figure('Position',[50,50,600,400]);
            angle=linspace(0,2*pi,100000);
            plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
            hold on
            for i =1:nsimul
                if t{i}(end) > 3.2e0
%                     plot(xT{i}.*cos(omega*t{i})+yT{i}.*sin(omega*t{i}), -xT{i}.*sin(omega*t{i})+yT{i}.*cos(omega*t{i}),'k+')
%                     hold on
%                     plot(xL{i}.*cos(omega*t{i})+yL{i}.*sin(omega*t{i}), -xL{i}.*sin(omega*t{i})+yL{i}.*cos(omega*t{i}),'r+')
%                     hold on
                    plot(xA{i}.*cos(omega*t{i})+yA{i}.*sin(omega*t{i}), -xA{i}.*sin(omega*t{i})+yA{i}.*cos(omega*t{i}))
                    hold on
                end
            end
            hold off
            xlabel("x' [m]")
            ylabel("y' [m]")
            axis equal 
            grid on
            set(gca,'fontsize',15);
            print(fig2,'figures/troisCorpsAdaptRho0_trajectoireRprime', '-depsc');
        end
    end 

    figure
    plot(theta, rentree,'+')
    xlabel('\theta_0 [rad]')
    ylabel('rentree (y/n)')
    grid on
    
    
    fig3=figure('Position',[50,50,600,400]);
    plot(theta,acc,'+')
    xlabel('\theta_0 [rad]')
    ylabel('max_{acc} [m/s^2]')
    grid on
    set(gca,'fontsize',15);
    print(fig3,'figures/troisCorpsAdaptRho0_acc', '-depsc');
    
    fig4=figure('Position',[50,50,600,400]);
    plot(theta,hmin-RT,'+')
    xlabel('\theta_0 [rad]')
    ylabel('min_{h} [m]')
    grid on
    set(gca,'fontsize',15);
    print(fig4,'figures/troisCorpsAdaptRho0_hmin', '-depsc');
    
    fig5=figure('Position',[50,50,600,400]);
    yyaxis left;
    plot(theta,hmin-RT,'+')
    hold on
    yyaxis right;
    plot(theta,acc,'+')
    hold on
    hold off
    yyaxis left
    xlabel('\theta_0 [rad]')
    ylabel('min_{h} [m]')
    yyaxis right
    ylabel('max_{acc} [m/s^2]')
    set(gca,'fontsize',13);
    grid on;
    print(fig5,'figures/troisCorpsAdaptRho0_hmin&acc', '-depsc');
    
    figure
    plot(t{76},Pt{76})
    xlabel('Temps [s]')
    ylabel('Puissance de traînée [J]')
    grid on

    figure
    plot(t{76},sqrt((xA{76}-xT{76}).^2+(yA{76}-yT{76}).^2))
    xlabel('Temps [s]')
    ylabel('Distance à la Terre [m]')
    grid on
    
    figure
    plot(t{76},vAbs{76},t{76},vRel{76})
    xlabel('Temps [s]')
    ylabel('Vitesse [m/s]')
    grid on
    
    
elseif strcmp(paraName, 'condIn')
    if trajectoire
        figure
        angle=linspace(0,2*pi,100000);
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


