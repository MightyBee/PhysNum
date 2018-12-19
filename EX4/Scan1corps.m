% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=1.2;
tFin=2*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             5000      % dt
             1e-6      % precision
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
variables = [0               0               r0              % x0
             0               384748000       0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             mT              7.3477e0        5809            % m
             RT              1737000         1.95            % R
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);


%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable

nsimul = 200; % Nombre de simulations à faire

% theta     = linspace(0,2*pi,nsimul);
% theta1    = [2.94 2.99; 3.3 3.35];
theta1    = [2.9505 2.952; 3.331 3.3325];
theta     = create_linspace(theta1,nsimul);
dt        = logspace(1.5,0,nsimul); % Valeurs du parametre a scanner
precision = logspace(-1,-8,nsimul); % Valeurs du parametre a scanner
% precision = (logspace(4,32,nsimul)).^(-1/4);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paraName='theta'; % Nom du parametre a scanner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


trajectoire=true;


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
    nsteps = ones(1,nsimul);
    if rho0==0
        hmin = zeros(1,nsimul);
        vmax = zeros(1,nsimul);
        dt= cell(2,3);
        r = cell(1,3);
    else
        maxAcc=zeros(1,nsimul);
        maxPT=zeros(1,nsimul);
    end
elseif strcmp(paraName, 'both')
    hmin   = zeros(2,nsimul/2);
    nsteps =  ones(2,nsimul/2);
elseif strcmp(paraName, 'theta')
%     error = zeros(1,nsimul);
%     maxAcc=zeros(1,nsimul);
%     maxPT=zeros(1,nsimul);
%     hmin = zeros(1,nsimul);
    t  = cell(1,nsimul);
    acc= zeros(1,nsimul);
    Pt = zeros(1,nsimul);
    xT = cell(1,nsimul);
    yT = cell(1,nsimul);
    xA = cell(1,nsimul);
    yA = cell(1,nsimul);
    hmin= zeros(1,nsimul);
    lmin= zeros(1,nsimul);
    vAbs= cell(1,nsimul);
    vRel= cell(1,nsimul);
    rentree=zeros(1,nsimul);
end

% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2  


for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
    if strcmp(paraName, 'dt') || strcmp(paraName, 'precision' )
        t = data(:,1);
        xA = data(:,17);
        yA = data(:,18);
        nsteps(i)=size(t,1)-1;
        if rho0==0
            xT = data(:,5);
            yT = data(:,6);
            vxT= data(:,8);
            vyT= data(:,9);
            vxA= data(:,20);
            vyA= data(:,21);
            hmin(i)=inter_min(t,sqrt((xA-xT).^2+(yA-yT).^2),3);
            vmax(i)=inter_max(t,sqrt((vxA-vxT).^2+(vyA-vyT).^2),3);
            if strcmp(paraName, 'precision')
                dt{1,ceil(i/nsimul*3)}=t(1:end-2);
                dt{2,ceil(i/nsimul*3)}=t(2:end-1)-t(1:end-2);
                nsteps(i)=2*nsteps(i);
                r{ceil(i/nsimul*3)}=sqrt((xA(1:end-2)-xT(1:end-2)).^2+(yA(1:end-2)-yT(1:end-2)).^2);
            end
        else
            maxAcc(i)= inter_max(t,data(:,2),3);
            maxPT(i) = inter_max(t,abs(data(:,4)),3);
            nsteps(i)=2*nsteps(i);
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
    elseif strcmp(paraName, 'theta')
%         t = data(:,1);
%         acc = data(:,2);
%         xT = data(:,5);
%         yT = data(:,6);
%         xA = data(:,17);
%         yA = data(:,18);
%         a=acc(1);
%         hmin(i);
%         for l=1:size(acc,1)
%             if acc(l)>a
%                 a=acc(l);
%             end
%         end
%         maxAcc(i)= a;
%         hmin(i)= sqrt((xA(1)-xT(1))^2+(yA(1)-yT(1))^2);
%         for j=1:size(t,1)
%             if sqrt((xA(j)-xT(j))^2+(yA(j)-yT(j))^2) < hmin(i)
%             hmin(i)=sqrt((xA(j)-xT(j))^2+(yA(j)-yT(j))^2);
%             end
%         end 
        t{i}  = data(:,1);
        acc(i)= max(data(:,2));
        Pt(i) = max(abs(data(:,4)));
        xT{i} = data(:,5);
        yT{i} = data(:,6);
        xA{i} = data(:,17);
        yA{i} = data(:,18);
        hmin(i)=min(sqrt((xA{i}-xT{i}).^2+(yA{i}-yT{i}).^2));
        vAbs{i}= sqrt(data(:,20).^2+data(:,21).^2);
        vRel{i}= sqrt((data(:,20)-data(:,8)).^2+(data(:,21)-data(:,9)).^2);
        if hmin(i)<=RT
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
            Pt(i) = inter_max(t{i},abs(data(1:iFin,4)),3);
%             Pt(i) = max(data(1:iFin,4));
            xT{i} = xT{i}(1:iFin);
            yT{i} = yT{i}(1:iFin);
            xA{i} = xA{i}(1:iFin);
            yA{i} = yA{i}(1:iFin);
            vAbs{i}= vAbs{i}(1:iFin);
            vRel{i}= vRel{i}(1:iFin);
        end
        if not (hmin(i)<=RT && rentree(i)==1)
            t{i}  = nan;
            acc(i)= nan;
            Pt(i) = nan;
            xT{i} = nan;
            yT{i} = nan;
            xA{i} = nan;
            yA{i} = nan;
            vAbs{i}= nan;
            vRel{i}= nan;
        end
    end
end


%% Figures %%
%%%%%%%%%%%%%
%

if strcmp(paraName, 'dt') || strcmp(paraName, 'precision')
    fig1=figure('Position',[50,50,700,420]);
    plot(RT*cos(linspace(0,2*pi,100000)),RT*sin(linspace(0,2*pi,100000)),'r')
    hold on
    plot(xA,yA)
    hold off
    ylim([-16e7 4e7])
    xlabel('x [m]')
    ylabel('y [m]')
    set(gca,'fontsize',15);
    axis equal
    grid on
    lgd=legend('Surface terrestre',"Trajectoire d'Appolo");
    set(lgd,'fontsize',14,'Location','southwest');
    rectangle('Position',[-6.3e6 1.8e6 0.5e6 1e6]);
    axes('position',[.60 .24 .25 .50])
    box on % put box around new pair of axes
    plot(RT*cos(linspace(0,2*pi,100000)),RT*sin(linspace(0,2*pi,100000)),'r',xA,yA) % plot on new axes
    xlim([-6.5e6 -5.8e6])
    ylim([0.95e6 1.9e6])
    axis equal
%     xlim([-6.3e6 -5.8e6])
%     ylim([1.8e6 2.4e6])
    grid on
    if rho0==0
        fig2=figure('Position',[50,50,600,400]);
        loglog(nsteps,abs(hmin-h-RT),'+',nsteps,abs(hmin(1)-h-RT)*nsteps(1)^4*nsteps.^(-4),'--')
        grid on
        xlabel('N_{steps}')
        ylabel('Erreur sur h_{min} [m]')
        set(gca,'fontsize',15);
        lgd=legend('Runge-Kutta 4','\propto 1/N^4');
        set(lgd,'fontsize',14,'Location','southwest');

        fig3=figure('Position',[50,50,600,400]);
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

        elseif strcmp(paraName, 'precision')

            fig4=figure('Position',[50,50,600,400]);
            plot(dt{1,1},dt{2,1},dt{1,2},dt{2,2},dt{1,3},dt{2,3})
            grid on
            xlabel('t [s]')
            ylabel('\Deltat [s]')
            set(gca,'fontsize',15);
            lgd=legend(sprintf('N_{steps}=%d',size(dt{1,1},1)+1),sprintf('N_{steps}=%d',size(dt{1,2},1)+1),sprintf('N_{steps}=%d',size(dt{1,3},1)+1));
            set(lgd,'fontsize',14,'Location','northeast');

            fig5=figure('Position',[50,50,600,400]);
            plot(r{1},dt{2,1},r{2},dt{2,2},r{3},dt{2,3})
            grid on
            xlabel('Distance Terre-Appolo [m]')
            ylabel('\Deltat [s]')
            set(gca,'fontsize',15);
            lgd=legend(sprintf('N_{steps}=%d',size(dt{1,1},1)+1),sprintf('N_{steps}=%d',size(dt{1,2},1)+1),sprintf('N_{steps}=%d',size(dt{1,3},1)+1));
            set(lgd,'fontsize',14,'Location','northwest');
            
            figure
            semilogx(precision,nsteps)

            print(fig1,'figures/unCorpsAdapt_trajectoire', '-depsc');
            print(fig2,'figures/unCorpsAdapt_convH', '-depsc');
            print(fig3,'figures/unCorpsAdapt_convV', '-depsc');
            print(fig4,'figures/unCorpsAdapt_dtT', '-depsc');
            print(fig5,'figures/unCorpsAdapt_dtR', '-depsc');
        end
    else
        fig2=figure('Position',[50,50,600,400]);
        plot(nsteps.^(-4) , maxAcc,'+')
        grid on
        xlabel('1/N_{steps}^4')
        ylabel('Acc_{max} [m/s^{2}]')
        set(gca,'fontsize',15);

        fig3=figure('Position',[50,50,600,400]);
        plot(nsteps.^(-4), maxPT,'+')
        grid on
        xlabel('1/N_{steps}^4')
        ylabel('P^{t}_{max} [W]')
        set(gca,'fontsize',15);
        
        print(fig1,'figures/unCorpsAdaptRho0_trajectoire', '-depsc');
        print(fig2,'figures/unCorpsAdaptRho0_convA', '-depsc');
        print(fig3,'figures/unCorpsAdaptRho0_convP', '-depsc');
    end

elseif strcmp(paraName , 'theta')

%     MAXA=[];
%     THETA=[];
%     for i=1:nsimul
%         r=0;
%         z=0;
%         if (hmin(i) <= RT)
%             for j=1:size(t,1);
%                 if sqrt((xA(j)-xT(j))^2+(yA(j)-yT(j))^2)==hmin(i);
%                     r=j;
%                 end
%             end
%           for l=1:(r-1)
%               if (xA(l+1)-xT(l+1))^2+(yA(l+1)-yT(j))^2>(xA(l)-xT(l))^2+(yA(l)-yT(l))^2;
%                   z=1;
%               end
%           end
%           if z==0;
%               MAXA=[MAXA maxAcc(i)] ;
%               THETA=[THETA theta(i)];
%           end
%         end
%     end
% 
%     figure
%     plot(THETA, MAXA ,'+')
%     grid on
%     
%     figure
%     plot(theta, hmin)
%     grid on
    if trajectoire
        fig1=figure('Position',[50,50,600,400]);
        angle=linspace(0,2*pi,100000);
        plot(RT*cos(angle),RT*sin(angle),'r')
        hold on
        for i =1:nsimul
            plot(xA{i},yA{i})
            hold on
        end
        hold off
        xlabel('x [m]')
        ylabel('y [m]')
        axis equal
        grid on
        set(gca,'fontsize',15);
        print(fig1,'figures/unCorpsAdaptRho0_theta_trajectoire', '-depsc');
    end

    figure
    plot(theta, rentree,'+')
    xlabel('\theta_0 [rad]')
    ylabel('rentree (y/n)')
    grid on
    
    fig2=figure('Position',[50,50,600,400]);
    plot(theta,hmin-RT,'p')
    xlabel('\theta_0 [rad]')
    ylabel('min(h) [m]')
    grid on    
    set(gca,'fontsize',15);
    print(fig2,'figures/unCorpsAdaptRho0_theta_hmin', '-depsc');
    
    fig3=figure('Position',[50,50,600,400]);
    plot(theta,acc,'+')
    xlabel('\theta_0 [rad]')
    ylabel('max(acc) [m/s^2]')
    grid on
    set(gca,'fontsize',15);
    print(fig3,'figures/unCorpsAdaptRho0_theta_acc', '-depsc');
    
    fig4=figure('Position',[50,50,600,400]);
    plot(theta,Pt,'+')
    xlabel('\theta_0 [rad]')
    ylabel('max(P_t) [W]')
    grid on
    set(gca,'fontsize',15);
    print(fig4,'figures/unCorpsAdaptRho0_theta_Pt', '-depsc');
    





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
