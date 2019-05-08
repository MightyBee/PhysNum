% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 
input = 'configuration.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nbEtude=10;
precision=logspace(-3,-5,nbEtude);

deltaMP=[50 150];

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

    
thetaFinal=zeros(1,nbEtude);
nstepsFinal=zeros(1,nbEtude);

for l = 1:nbEtude
    nsimul = 2;
    
    deltaM=deltaMP(1);
    deltaP=deltaMP(2);
    visee=0.5;
    erreur=1e-2;
    
    paramstr = {"delta"; "x0"};
    param = [deltaM deltaP]; % Valeurs du parametre a scanner
    
    output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
    for i = 1:nsimul
        % Sring des paramètres à varier
        parameter = '';
        for k=1:size(paramstr,1)
            parameter=[parameter sprintf('%s=%.15g ', paramstr{k}, param(k,i))];
        end
        parameter=strip(parameter);
        % Nom du fichier de sortie
        output{i} = [dossier strrep(parameter, ' ', '_')];
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s output=%s', repertoire, executable, input, parameter, output{i});
        disp(cmd)
        system(cmd);
    end
    
    deltata0=(deltaM+deltaP)/2;
    
    err=[0 0];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:2
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_g{i}=data(:,2);
        prob_d{i}=data(:,3);
        prob_dmax(i)=max(data(:,3));
        r=sqrt((xA-xT).^2+(yA-xT).^2);
        [value, indice]=min(r);
        if indice+3<=size(t,1)
            err(i)=inter_min(t,r,3)-RT;
        end
    end
    
    hMinM=err(1);
    hMinP=err(2);
    
    %% Dichotomie %%
    
    while max(hMinM,hMinP)>visee && min(hMinM,hMinP)<visee
        
        paramstr = {"vx0"; "vy0"};
        param = [v0A*cos(deltata0); v0A*sin(deltata0)+omega*rA]; % Valeurs du parametre a scanner
        input=sprintf('configuration%d.in', configfileNb);
        parameter = "";
        for i=1:size(paramstr,1)
            parameter=parameter+sprintf('%s=%.20g ', paramstr{i}, param(i));
        end
        parameter=strip(parameter);
        % Nom du fichier de sortie
        output = "simulations/"+strrep(parameter, ' ', '_')+sprintf('_precision=%.15g.out',precision(l));
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %d %s configuration0.in 3 tFin=%.15g precision=%.15g output=%s', repertoire, executable, input, size(param,1), parameter, tFin, precision(l), output);
        disp(cmd)
        system(cmd);
        
        % 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
        % t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3
        
        data = load(output); % Chargement du fichier de sortie de la i-ieme simulation
        t=data(:,1);
        xT = data(:,5);
        yT = data(:,6);
        xL = data(:,11);
        yL = data(:,12);
        xA = data(:,17);
        yA = data(:,18);
        [value, indice]=min(sqrt((xA-xT).^2+(yA-yT).^2))
        if indice+3>size(t,1)
            hMin0=0;
        else
            hMin0=inter_min(t,sqrt((xA-xT).^2+(yA-yT).^2),3)-RT
        end
        
        theta=[theta deltata0];
        err=[err hMin0];
        nsteps=size(t,1)-1;
        
        disp([deltaM deltata0 deltaP]);
        disp([hMinM hMin0 hMinP]);
        
        if abs(hMin0-visee) < erreur
            disp('OKAY')
            break;
        else
            if sign(hMinP-hMinM)*(hMin0-visee) > erreur
                if hMin0==0
                    deltaP=deltata0;
                    hMinP=hMin0;
                    deltata0=(9*deltaM+deltaP)/10;
                elseif hMinP==0
                    deltaP=deltata0;
                    hMinP=hMin0;
                    deltata0=deltaM-(hMinM-visee)*(deltaP-deltaM)/(hMinP-hMinM);
                else
                    theta1=(hMinP*deltata0-hMin0*deltaP)/(hMinP-hMin0);
                    deltaP=deltata0;
                    hMinP=hMin0;
                    deltata0=theta1+visee*(deltata0-theta1)/hMin0;
                end
            elseif sign(hMinP-hMinM)*(hMin0-visee) < -erreur
                if hMin0==0
                    deltaM=deltata0;
                    hMinM=hMin0;
                    deltata0=(deltaM+9*deltaP)/10;
                elseif hMinM==0
                    deltaM=deltata0;
                    hMinM=hMin0;
                    deltata0=deltaM-(hMinM-visee)*(deltaP-deltaM)/(hMinP-hMinM);
                else
                    theta1=(hMinM*deltata0-hMin0*deltaM)/(hMinM-hMin0);
                    deltaM=deltata0;
                    hMinM=hMin0;
                    deltata0=theta1+visee*(deltata0-theta1)/hMin0;
                end
            else
                warning('ERRRRREEEEUUUUURRRRRR');
                break;
            end
        end
    end
    disp(num2str(deltata0,15))
    disp(hMin0)
    disp(nsteps)
    thetaFinal(l)=deltata0;
    nstepsFinal(l)=nsteps;
end

%% Analyse %%
%%%%%%%%%%%%%
for nb=1:nTheta
    [p,~,mu]=polyfit(1./nstepsFinal.^4,thetaFinal,1);
    xFit=linspace(0,max(1./nstepsFinal.^4),100000);
    yFit=polyval(p,xFit,[],mu);
end
%% Figures %%
%%%%%%%%%%%%%

fig=[];

for nb=1:nTheta
    fig=[fig figure('Position',[50,50,600,400])];
    plot(1./nstepsFinal.^4,thetaFinal,'b+')
    hold on
    plot(xFit,yFit,'r--',xFit(1),yFit(1),'rp', 'MarkerSize',2)
    hold off
    grid on
    xlabel('1/N_{steps}^4')
    ylabel('\theta_0 [rad]')
    set(gca,'fontsize',15);
    lgd=legend('Simulations','Régression linéaire');
    if nb==1 || nb==3
        set(lgd,'fontsize',14,'Location','southeast');
    else
        set(lgd,'fontsize',14,'Location','northeast');
    end
    print(fig(nb),"figures/troisCorps_convTheta0_"+num2str(nb), '-depsc');
end

fig2=figure('Position',[50,50,600,400]);
for nb=1:nTheta
    plot(xT,yT,'r',xL,yL,'k',xA,yA);
    hold on
end
hold off
grid on
axis equal
xlabel('x [m]')
ylabel('y [m]')
xlim([-0.5e8 5.2e8]);
ylim([-1.5e8, 2.5e8]);
set(gca,'fontsize',15);
% lgd=legend('Trajectoire de la Terre','Trajectoire de la Lune',"Trajectoire d'Appolo");
% set(lgd,'fontsize',14,'Location','northwest');
print(fig2,'figures/troisCorps_trajTheta0', '-depsc');

fig3=figure('Position',[50,50,600,400]);
plot(rT+RT*cos(linspace(0,2*pi,10000)),RT*sin(linspace(0,2*pi,10000)),'r')
hold on
plot(rL+RL*cos(linspace(0,2*pi,10000)),RL*sin(linspace(0,2*pi,10000)),'k')
hold on
for nb =1:nTheta
    plot(xA.*cos(omega*t)+yA.*sin(omega*t), -xA.*sin(omega*t)+yA.*cos(omega*t))
    hold on
end 
hold off
grid on
axis equal
xlabel("x' [m]")
ylabel("y' [m]")
set(gca,'fontsize',15);
% lgd=legend('Surface terrestre','Surface lunaire',"Trajectoire d'Appolo");
% set(lgd,'fontsize',14,'Location','southeast');
print(fig3,'figures/troisCorps_trajTheta0Prime', '-depsc');

for nb=1:nTheta
    fprintf('vx0=%.15g \n vy0=%.15g \n',v0A*cos(yFit(1)),v0A*sin(yFit(1)));
end
