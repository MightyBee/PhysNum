% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%

G=6.674e-11;
rho0=0;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             2*24*3600 % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             50        % dt
             1e-5      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             1      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);

rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

alpha = 0.3; %0.2138; %pi-asin(vMax_th*(h+RT)/(v0A*r0));
v0A=1200;
vx0A=v0A*cos(alpha);
vy0A=v0A*sin(alpha)+omega*rA;

vy0T=omega*rT;
vy0L=omega*rL;

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
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);



%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 

nbEtude=10;
precision=logspace(-3,-5,nbEtude);


tFin=   [2.5;     3;       6.5;   14;       14.5     ]*24*3600;
thetaMP=[3.6 3.8; 4.2 4.4; 5 5.2; 5.6 5.63; 0.28 0.29];
nTheta = size(thetaMP,1);


nstepsFinal=cell(1,nTheta);
thetaFinal=cell(1,nTheta);
xFit=cell(1,nTheta);
yFit=cell(1,nTheta);
xT=cell(1,nTheta);
yT=cell(1,nTheta);
xL=cell(1,nTheta);
yL=cell(1,nTheta);
xA=cell(1,nTheta);
yA=cell(1,nTheta);
t=cell(1,nTheta);


for nb=1:nTheta
    
    thetaFinal{nb}=zeros(1,nbEtude);
    nstepsFinal{nb}=zeros(1,nbEtude);

    for l = 1:nbEtude
        nsimul = 2;

        thetaM=thetaMP(nb,1);
        thetaP=thetaMP(nb,2);
        visee=10000;
        erreur=1;

        paramstr = {"vx0"; "vy0"};
        theta = [thetaM thetaP];
        param = [v0A*cos(theta); v0A*sin(theta)+omega*rA]; % Valeurs du parametre a scanner
        configfileNb=3;

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
            output{1,k} = "simulations/"+strrep(parameter, ' ', '_')+sprintf('_nbTheta=%d_precision=%.15g.out',nb,precision(l));
            % Execution du programme en lui envoyant la valeur a scanner en argument
            cmd = sprintf('%s%s %s %d %s configuration0.in 3 tFin=%.15g precision=%.15g output=%s', repertoire, executable, input, size(param,1), parameter, tFin(nb), precision(l), output{1,k});
            disp(cmd)
            system(cmd);
        end


        theta0=(thetaM+thetaP)/2;

        % 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
        % t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3

        hmin=[0 0];

%         figure
        for i=1:2
            data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
            t{nb}  = data(:,1);
            xT{nb} = data(:,5);
            yT{nb} = data(:,6);
            xA{nb} = data(:,17);
            yA{nb} = data(:,18);
            r=sqrt((xA{nb}-xT{nb}).^2+(yA{nb}-xT{nb}).^2);
            [value, indice]=min(r);
            if indice+3<=size(t{nb},1)
                hmin(i)=inter_min(t{nb},r,3)-RT;
            end
%             plot(xA{nb},yA{nb})
%             hold on
        end
%         hold off
%         legend(num2str(hmin(1)),num2str(hmin(2)))

        hMinM=hmin(1);
        hMinP=hmin(2);

        %% Dichotomie %%

        while max(hMinM,hMinP)>visee && min(hMinM,hMinP)<visee

            paramstr = {"vx0"; "vy0"};
            param = [v0A*cos(theta0); v0A*sin(theta0)+omega*rA]; % Valeurs du parametre a scanner
            input=sprintf('configuration%d.in', configfileNb);
            parameter = "";
            for i=1:size(paramstr,1)
                parameter=parameter+sprintf('%s=%.20g ', paramstr{i}, param(i));
            end
            parameter=strip(parameter);
            % Nom du fichier de sortie
            output = "simulations/"+strrep(parameter, ' ', '_')+sprintf('_nbTheta=%d_precision=%.15g.out',nb,precision(l));
            % Execution du programme en lui envoyant la valeur a scanner en argument
            cmd = sprintf('%s%s %s %d %s configuration0.in 3 tFin=%.15g precision=%.15g output=%s', repertoire, executable, input, size(param,1), parameter, tFin(nb), precision(l), output);
            disp(cmd)
            system(cmd);

            % 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
            % t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3

            data = load(output); % Chargement du fichier de sortie de la i-ieme simulation
            t{nb}=data(:,1);
            xT{nb} = data(:,5);
            yT{nb} = data(:,6);
            xL{nb} = data(:,11);
            yL{nb} = data(:,12);
            xA{nb} = data(:,17);
            yA{nb} = data(:,18);
            [value, indice]=min(sqrt((xA{nb}-xT{nb}).^2+(yA{nb}-yT{nb}).^2));
            if indice+3>size(t,1)
                hMin0=0;
            else
                hMin0=inter_min(t{nb},sqrt((xA{nb}-xT{nb}).^2+(yA{nb}-yT{nb}).^2),3)-RT;
            end

            theta=[theta theta0];
            hmin=[hmin hMin0];
            nsteps=size(t{nb},1)-1;

            disp([thetaM theta0 thetaP]);
            disp([hMinM hMin0 hMinP]);

            if abs(hMin0-visee) < erreur
                disp('OKAY')
                break;
            else
                if sign(hMinP-hMinM)*(hMin0-visee) > erreur
                    if hMin0==0
                        thetaP=theta0;
                        hMinP=hMin0;
                        theta0=(9*thetaM+thetaP)/10;
                    elseif hMinP==0
                        thetaP=theta0;
                        hMinP=hMin0;
                        theta0=thetaM-(hMinM-visee)*(thetaP-thetaM)/(hMinP-hMinM);
                    else
                        theta1=(hMinP*theta0-hMin0*thetaP)/(hMinP-hMin0);
                        thetaP=theta0;
                        hMinP=hMin0;
                        theta0=theta1+visee*(theta0-theta1)/hMin0;
                    end
                elseif sign(hMinP-hMinM)*(hMin0-visee) < -erreur
                    if hMin0==0
                        thetaM=theta0;
                        hMinM=hMin0;
                        theta0=(thetaM+9*thetaP)/10;
                    elseif hMinM==0
                        thetaM=theta0;
                        hMinM=hMin0;
                        theta0=thetaM-(hMinM-visee)*(thetaP-thetaM)/(hMinP-hMinM);
                    else
                        theta1=(hMinM*theta0-hMin0*thetaM)/(hMinM-hMin0);
                        thetaM=theta0;
                        hMinM=hMin0;
                        theta0=theta1+visee*(theta0-theta1)/hMin0;
                    end
                else 
                    warning('ERRRRREEEEUUUUURRRRRR');
                    break;
                end
            end
        end
        disp(num2str(theta0,15))
        disp(hMin0)
        disp(nsteps)
        thetaFinal{nb}(l)=theta0;
        nstepsFinal{nb}(l)=nsteps;
    end
end

%% Analyse %%
%%%%%%%%%%%%%
for nb=1:nTheta
    [p,~,mu]=polyfit(1./nstepsFinal{nb}.^4,thetaFinal{nb},1);
    xFit{nb}=linspace(0,max(1./nstepsFinal{nb}.^4),100000);
    yFit{nb}=polyval(p,xFit{nb},[],mu);
end
%% Figures %%
%%%%%%%%%%%%%

fig=[];

for nb=1:nTheta
    fig=[fig figure('Position',[50,50,500,400])];
    plot(1./nstepsFinal{nb}.^4,thetaFinal{nb},'b+')
    hold on
    plot(xFit{nb},yFit{nb},'r-',xFit{nb}(1),yFit{nb}(1),'rp', 'MarkerSize',2)
    hold off
    grid on
    xlabel('1/N_{steps}^4')
    ylabel('\theta_0 [rad]')
    set(gca,'fontsize',15);
    lgd=legend('Simulations','Régression linéaire');
    set(lgd,'fontsize',14,'Location','northwest');
    print(fig(nb),'figures/troisCorps_convTheta0', '-depsc');
end 

fig2=figure('Position',[50,50,600,400]);
for nb=1:nTheta
    plot(xT{nb},yT{nb},'r',xL{nb},yL{nb},'k',xA{nb},yA{nb});
    hold on
end
hold off
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
ylim([-0.8e8 2.2e8]);
set(gca,'fontsize',15);
lgd=legend('Trajectoire de la Terre','Trajectoire de la Lune',"Trajectoire d'Appolo");
set(lgd,'fontsize',14,'Location','northwest');
print(fig2,'figures/troisCorps_trajTheta0', '-depsc');

fig3=figure('Position',[50,50,600,400]);
plot(rT+RT*cos(linspace(0,2*pi,10000)),RT*sin(linspace(0,2*pi,10000)),'r')
hold on
plot(rL+RL*cos(linspace(0,2*pi,10000)),RL*sin(linspace(0,2*pi,10000)),'k')
hold on
for nb =1:nTheta
    plot(xA{nb}.*cos(omega*t{nb})+yA{nb}.*sin(omega*t{nb}), -xA{nb}.*sin(omega*t{nb})+yA{nb}.*cos(omega*t{nb}))
    hold on
end 
hold off
grid on
axis equal
xlabel("x' [m]")
ylabel("y' [m]")
ylim([-2.5e8 0.5e8]);
set(gca,'fontsize',15);
lgd=legend('Surface terrestre','Surface lunaire',"Trajectoire d'Appolo");
set(lgd,'fontsize',14,'Location','southeast');
print(fig3,'figures/troisCorps_trajTheta0Prime', '-depsc');

for nb=1:nTheta
    fprintf('vx0=%.15g \n vy0=%.15g',v0A*cos(yFit{nb}(1)),v0A*sin(yFit{nb}(1)));
end
