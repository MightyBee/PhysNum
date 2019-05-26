%% Initialisation %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'Exercice8'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'classique.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nbEtude=1;
dt=5*logspace(0,-1,nbEtude);
dt=1;

thetaMP=[11 11.1];



thetaFinal=zeros(1,nbEtude);
nstepsFinal=zeros(1,nbEtude);

for l = 1:nbEtude
    nsimul = 2;
    
    thetaM=thetaMP(1);
    thetaP=thetaMP(2);
    visee=0.5;
    erreur=1e-6;
    
    paramstr = {"n"; "dt"};
    theta = [thetaM thetaP];
    param = [theta; dt(l)*ones(size(theta))]; % Valeurs du parametre a scanner
    
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
        cmd = sprintf('%s%s %s %s delta=64 x0=-64 output=%s', repertoire, executable, input, parameter, output{i});
        disp(cmd)
        system(cmd);
    end
    
    
    Pmax=[0 0];
    
    for i=1:2
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_d=data(:,3);
        Pmax(i)=max(prob_d(1:round(2*length(t)/5)));
    end
    
    PmaxM=Pmax(1)
    PmaxP=Pmax(2)
    
    theta0=thetaP-(PmaxP-visee)*(thetaP-thetaM)/(PmaxP-PmaxM);% (thetaM+thetaP)/2;
    
    %% Dichotomie %%
    
    while max(PmaxM,PmaxP)>visee && min(PmaxM,PmaxP)<visee
        
        paramstr = {"n", "dt"};
        param = [theta0 dt(l)]; % Valeurs du parametre a scanner
        parameter = '';
        for k=1:size(paramstr,1)
            parameter=[parameter sprintf('%s=%.15g ', paramstr{k}, param(k))];
        end
        parameter=strip(parameter);
        % Nom du fichier de sortie
        output = [dossier strrep(parameter, ' ', '_')];
        % Execution du programme en lui envoyant la valeur a scanner en argument
        cmd = sprintf('%s%s %s %s delta=64 x0=-64 output=%s', repertoire, executable, input, parameter, output);
        disp(cmd)
        system(cmd);
        
        
        data = load([output,'_obs.out']);
        t=data(:,1);
        prob_d=data(:,3);
        Pmax0=max(prob_d(1:round(2*length(t)/5)));
        
        theta=[theta theta0];
        Pmax=[Pmax Pmax0];
        
        disp([thetaM theta0 thetaP]);
        disp([PmaxM Pmax0 PmaxP]);
        
        if abs(Pmax0-visee) < erreur
            disp('OKAY')
            fprintf('erreur : %.4g\n\n', abs(Pmax0-visee))
            break;
        else
            if sign(PmaxP-PmaxM)*(Pmax0-visee) > erreur
                        theta1=(PmaxP*theta0-Pmax0*thetaP)/(PmaxP-Pmax0);
                    thetaP=theta0;
                    PmaxP=Pmax0;
                    theta0=theta1+visee*(theta0-theta1)/Pmax0;
                
            elseif sign(PmaxP-PmaxM)*(Pmax0-visee) < -erreur
                
                    theta1=(PmaxM*theta0-Pmax0*thetaM)/(PmaxM-Pmax0);
                    thetaM=theta0;
                    PmaxM=Pmax0;
                    theta0=theta1+visee*(theta0-theta1)/Pmax0;
                
            else
                warning('ERRRRREEEEUUUUURRRRRR');
                break;
            end
        end
    end
    disp(num2str(theta0,15))
    disp(Pmax0)
    thetaFinal(l)=theta0;
end

%%

fig1=figure('Position',[50,50,600,450]);
semilogy(theta, abs(Pmax-visee),'r+','MarkerSize',11)
set(gca,'FontSize',15)
xlabel('$n$','Interpreter','Latex','FontSize',22)
ylabel('$\left| P_{\rm trans}-0.5 \right|$','Interpreter','Latex','FontSize',22)
grid on, box off
print(fig1,'figures/n_dichotomie', '-depsc');
    
figure
semilogy(theta, abs(Pmax-visee),'+')
grid on
dd=[0 0:(length(Pmax)-2)];
figure
plot(dd,theta,'+')

% %% Analyse %%
% %%%%%%%%%%%%%
% [p,~,mu]=polyfit(1./nstepsFinal.^4,thetaFinal,1);
% xFit=linspace(0,max(1./nstepsFinal.^4),100000);
% yFit=polyval(p,xFit,[],mu);
% 
% %% Figures %%
% %%%%%%%%%%%%%
% 
% fig=[];
% 
% 
% fig=[fig figure('Position',[50,50,600,400])];
% plot(1./nstepsFinal.^4,thetaFinal,'b+')
% hold on
% plot(xFit,yFit,'r--',xFit(1),yFit(1),'rp', 'MarkerSize',2)
% hold off
% grid on
% xlabel('1/N_{steps}^4')
% ylabel('\theta_0 [rad]')
% set(gca,'fontsize',15);
% lgd=legend('Simulations','Régression linéaire');
% if nb==1 || nb==3
%     set(lgd,'fontsize',14,'Location','southeast');
% else
%     set(lgd,'fontsize',14,'Location','northeast');
% end
% print(fig,"figures/troisCorps_convTheta0", '-depsc');
% 
% fig2=figure('Position',[50,50,600,400]);
% plot(xT,yT,'r',xL,yL,'k',xA,yA);
% hold off
% grid on
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')
% xlim([-0.5e8 5.2e8]);
% ylim([-1.5e8, 2.5e8]);
% set(gca,'fontsize',15);
% % lgd=legend('Trajectoire de la Terre','Trajectoire de la Lune',"Trajectoire d'Appolo");
% % set(lgd,'fontsize',14,'Location','northwest');
% print(fig2,'figures/troisCorps_trajTheta0', '-depsc');
% 
% fig3=figure('Position',[50,50,600,400]);
% hold on
% plot(rT+RT*cos(linspace(0,2*pi,10000)),RT*sin(linspace(0,2*pi,10000)),'r')
% plot(rL+RL*cos(linspace(0,2*pi,10000)),RL*sin(linspace(0,2*pi,10000)),'k')
% plot(xA.*cos(omega*t)+yA.*sin(omega*t), -xA.*sin(omega*t)+yA.*cos(omega*t))
% hold off
% grid on
% axis equal
% xlabel("x' [m]")
% ylabel("y' [m]")
% set(gca,'fontsize',15);
% % lgd=legend('Surface terrestre','Surface lunaire',"Trajectoire d'Appolo");
% % set(lgd,'fontsize',14,'Location','southeast');
% print(fig3,'figures/troisCorps_trajTheta0Prime', '-depsc');
% 
% 
% fprintf('vx0=%.15g \n vy0=%.15g \n',v0A*cos(yFit(1)),v0A*sin(yFit(1)));
% 
