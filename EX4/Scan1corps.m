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
             3000      % dt
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
alpha = pi-asin(vMax_th*(h+RT)/(v0*r0))
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraName='precision'; % Nom du parametre a scanner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if strcmp(paraName,'theta')
    paramstr = {"vx0"; "vy0"};
    v=1200;
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
elseif strcmp(paraName,'both')
    paramstr = {"dt"; "precision"; "adaptatif"};
    param = [dt dt; precision precision; zeros(1,nsimul) ones(1,nsimul)];
    nsimul=2*nsimul;
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
    dt= cell(2,3);  
    r = cell(1,3);
elseif strcmp(paraName, 'both')
    hmin   = zeros(2,nsimul/2);
    nsteps =  ones(2,nsimul/2);
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
        vxT= data(:,8);
        vyT= data(:,9);
        xA = data(:,17);
        yA = data(:,18);
        vxA= data(:,20);
        vyA= data(:,21);
        nsteps(i)=size(t,1)-1;
        R = T1.Terre(8);
        hmin(i)=inter_min(t,sqrt((xA-xT).^2+(yA-yT).^2),3);
        vmax(i)=inter_max(t,sqrt((vxA-vxT).^2+(vyA-vyT).^2),3);
        if strcmp(paraName, 'precision')
            dt{1,ceil(i/nsimul*3)}=t(1:end-2);
            dt{2,ceil(i/nsimul*3)}=t(2:end-1)-t(1:end-2);
            nsteps(i)=2*nsteps(i);
            r{ceil(i/nsimul*3)}=sqrt((xA(1:end-2)-xT(1:end-2)).^2+(yA(1:end-2)-yT(1:end-2)).^2);
        end
    elseif strcmp(paraName,'both')
        t = data(:,1);
        xT = data(:,5);
        yT = data(:,6);
        xA = data(:,17);
        yA = data(:,18);
        if i>nsimul/2
            nsteps(2,i-nsimul/2)=2*(size(t,1)-1);
            hmin(2,i-nsimul/2)=inter_min(t,sqrt((xA-xT).^2+(yA-yT).^2),3);
        else
            nsteps(1,i)=size(t,1)-1;
            hmin(1,i)=inter_min(t,sqrt((xA-xT).^2+(yA-yT).^2),3);
        end
    elseif strcmp(paramstr, 'theta')
        t = data(:,1);
    end
end


%% Figures %%
%%%%%%%%%%%%%
%

if strcmp(paraName, 'dt') || strcmp(paraName, 'precision')
    fig1=figure('Position',[50,50,600,400])
    plot(RT*cos(linspace(0,2*pi,100000)),RT*sin(linspace(0,2*pi,100000)),'r')
    hold on
    plot(xA,yA)
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    set(gca,'fontsize',15);
    axis equal
    grid on
    lgd=legend('Surface terrestre',"Trajectoire d'Appolo");
    set(lgd,'fontsize',14,'Location','northwest');
 

    fig2=figure('Position',[50,50,600,400])
    loglog(nsteps,abs(hmin-h-RT),'+',nsteps,abs(hmin(1)-h-RT)*nsteps(1)^4*nsteps.^(-4),'--')
    grid on
    xlabel('N_{steps}')
    ylabel('Erreur sur h_{min} [m]')
    set(gca,'fontsize',15);
    lgd=legend('Runge-Kutta 4','\propto 1/N^4');
    set(lgd,'fontsize',14,'Location','southwest');

    fig3=figure('Position',[50,50,600,400])
    loglog(nsteps,abs(vmax-vMax_th),'+',nsteps,abs(vmax(1)-vMax_th)*nsteps(1)^4*nsteps.^(-4),'--')
    grid on
    xlabel('N_{steps}')
    ylabel('Erreur sur v_{max} [m/s]')
    set(gca,'fontsize',15);
    lgd=legend('Runge-Kutta 4','\propto 1/N^4');
    set(lgd,'fontsize',14,'Location','southwest');

    if strcmp(paraName, 'dt')
        print(fig1,'figures/unCorpsFixe_trajectoire', '-depsc');
        print(fig2,'figures/unCorpsFixe_convH', '-depsc');
        print(fig3,'figures/unCorpsFixe_convV', '-depsc');  
    else
        fig4=figure('Position',[50,50,600,400])
        plot(dt{1,1},dt{2,1},dt{1,2},dt{2,2},dt{1,3},dt{2,3})
        grid on
        xlabel('t [s]')
        ylabel('\Deltat [s]')
        set(gca,'fontsize',15);
        lgd=legend(sprintf('N_{steps}=%d',size(dt{1,1},1)+1),sprintf('N_{steps}=%d',size(dt{1,2},1)+1),sprintf('N_{steps}=%d',size(dt{1,3},1)+1));
        set(lgd,'fontsize',14,'Location','northeast');
        
        fig5=figure('Position',[50,50,600,400])
        plot(r{1},dt{2,1},r{2},dt{2,2},r{3},dt{2,3})
        grid on
        xlabel('Distance Terre-Appolo [m]')
        ylabel('\Deltat [s]')
        set(gca,'fontsize',15);
        lgd=legend(sprintf('N_{steps}=%d',size(dt{1,1},1)+1),sprintf('N_{steps}=%d',size(dt{1,2},1)+1),sprintf('N_{steps}=%d',size(dt{1,3},1)+1));
        set(lgd,'fontsize',14,'Location','northwest');
        
        print(fig1,'figures/unCorpsAdapt_trajectoire', '-depsc');
        print(fig2,'figures/unCorpsAdapt_convH', '-depsc');
        print(fig3,'figures/unCorpsAdapt_convV', '-depsc');
        print(fig4,'figures/unCorpsAdapt_dtT', '-depsc');
        print(fig5,'figures/unCorpsAdapt_dtR', '-depsc');    
    end
elseif strcmp(paraName,'both')
    fig1=figure('Position',[50,50,600,400]);
    disp(nsteps)
    disp(hmin)
    loglog(abs(hmin(1,:)-h-RT),nsteps(1,:),'+',abs(hmin(2,:)-h-RT),nsteps(2,:),'+')%,nsteps,abs(hmin(1)-h-RT)*nsteps(1)^4*nsteps.^(-4),'--')
    grid on
    xlabel('Erreur sur h_{min} [m]')
    ylabel('N_{steps}')
    set(gca,'fontsize',15);
    lgd=legend('Runge-Kutta 4 fixe','Runge-Kutta 4 adaptatif');
    set(lgd,'fontsize',14,'Location','southwest');
    print(fig1,'figures/unCorps_compFixAdapt', '-depsc');
end

