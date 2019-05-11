% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.


%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice8'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 11; % Nombre de simulations a faire
paraName='V0';

m=1;
omega=0.003;
P0=2*pi*14/400;
E0=P0^2/2; %mean(data(:,4));

if strcmp(paraName,'dt')
    dt=logspace(-0.2,1, nsimul);
    paramstr = {"dt"}; % Nom du parametre a scanner
    param = dt; % Valeurs du parametre a scanner
elseif strcmp(paraName,'Ninters')
    Ninters=logspace(2.1,3.1,nsimul);
    paramstr = {"Ninters"}; % Nom du parametre a scanner
    param = [Ninters]; % Valeurs du parametre a scanner
elseif strcmp(paraName,'delta')
    delta=linspace(0,150,nsimul);
    paramstr = {"delta"; "x0"}; % Nom du parametre a scanner
    param = [delta; -delta]; % Valeurs du parametre a scanner
elseif strcmp(paraName,'V0')
    V0=[0.1 1 3.1]*E0;
    nsimul=3;
    delta=sqrt(2*V0/(m*omega^2));
    paramstr = {"delta"; "x0"}; % Nom du parametre a scanner
    param = [delta; -delta]; % Valeurs du parametre a scanner
elseif strcmp(paraName,'n')
    n0=9;
    n1=13;
    nsimul=n1-n0+1;
    n=round(linspace(n0,n1,nsimul));
    delta = 64*ones(1,nsimul);
    paramstr = {"n"; "delta"; "x0"}; % Nom du parametre a scanner
    param = [n; delta; -delta]; % Valeurs du parametre a scanner
end


%% Simulations %%
%%%%%%%%%%%%%%%%%

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

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

if(strcmp(paraName,'dt'))
    xmoy = zeros(1,nsimul);
    pmoy = zeros(1,nsimul);
    dxmoy = zeros(1,nsimul);
    dpmoy = zeros(1,nsimul);
    errx = zeros(1,nsimul);
    errp = zeros(1,nsimul);
elseif(strcmp(paraName,'Ninters'))
    xmoy = zeros(1,nsimul);
    pmoy = zeros(1,nsimul);
    dxmoy = zeros(1,nsimul);
    dpmoy = zeros(1,nsimul);
    errx = zeros(1,nsimul);
    errp = zeros(1,nsimul);
elseif(strcmp(paraName,'delta'))
    prob_g = cell(1,nsimul);
    prob_d = cell(1,nsimul);
    prob_dmax = zeros(1,nsimul);
    V0=0.5*omega^2*delta.^2;
elseif(strcmp(paraName,'V0'))
    prob_g = cell(1,nsimul);
    prob_d = cell(1,nsimul);
    prob_dmax = zeros(1,nsimul);
    psi2 = cell(1,nsimul);
elseif(strcmp(paraName,'n'))
    prob_d = cell(1,nsimul);
end

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    if(strcmp(paraName,'dt'))
        data = load([output{i} '_obs.out']);
        t=data(:,1);
        x=data(:,5);
        xmoy(i) = x(end);
        x2moy   = data(end,6);
        dxmoy(i) = sqrt(x2moy-xmoy(i)^2);
        p=data(:,7);
        pmoy(i) = p(end);
        p2moy   = data(end,8);
        dpmoy(i) = sqrt(p2moy-pmoy(i)^2);
        A=sqrt(2*E0)/omega; 
%         A=max(x);
        errx(i)=max(abs(x/A-sin(omega*t)));
        errp(i)=max(abs(p/(A*omega)-cos(omega*t)));
    elseif(strcmp(paraName,'Ninters'))
        data = load([output{i} '_obs.out']);
        t=data(:,1);
        x=data(:,5);
        xmoy(i) = x(end);
        x2moy   = data(end,6);
        dxmoy(i) = sqrt(x2moy-xmoy(i)^2);
        p=data(:,7);
        pmoy(i) = p(end);
        p2moy   = data(end,8);
        dpmoy(i) = sqrt(p2moy-pmoy(i)^2);
        A=sqrt(2*E0)/omega; 
%         A=max(x);
        errx(i)=max(abs(x/A-sin(omega*t)));
        errp(i)=max(abs(p/(A*omega)-cos(omega*t)));
    elseif(strcmp(paraName,'delta'))
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_g{i}=data(:,2);
        prob_d{i}=data(:,3);
        prob_dmax(i)=max(data(:,3));
        E0=mean(data(:,4));
    elseif(strcmp(paraName,'V0'))
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_g{i}=data(:,2);
        prob_d{i}=data(:,3);
        prob_dmax(i)=max(data(:,3));
        data = load([output{i},'_pot.out']);
        x = data(:,1);
        psi2{i} = load([output{i},'_psi2.out']);
    elseif(strcmp(paraName,'n'))
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_d{i}=data(:,3);
    end
end


%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paraName,'dt'))
    fig1=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,xmoy,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\langle x \rangle \ \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig1,'figures/conv0_dt_x', '-depsc');

    fig2=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,dxmoy,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\langle \Delta x \rangle \ \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig2,'figures/conv0_dt_dx', '-depsc');
    
    fig3=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,pmoy,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\langle p \rangle \ \rm [kg \ m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig3,'figures/conv0_dt_p', '-depsc');

    fig4=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,dpmoy,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\langle \Delta p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig4,'figures/conv0_dt_dp', '-depsc');

    fig5=figure('Position',[50,50,600,450]);
    h=plot(dt.^2, errx,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0}|} \quad \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig5,'figures/conv0_dt_x_th', '-depsc');

    fig6=figure('Position',[50,50,600,450]);
    h=plot(dt.^2, errp,'k+');
    xlabel('$(\Delta t)^2 \ \rm [s^2]$','Interpreter','Latex')
    ylabel('$\max_t{|\frac{\langle p \rangle - p_{\rm th}}{A_0 \omega}|} \quad \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig6,'figures/conv0_dt_p_th', '-depsc');
    
   
elseif(strcmp(paraName,'Ninters'))
    fig1=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,xmoy,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\langle x \rangle \ \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig1,'figures/conv0_dx_x', '-depsc');

    fig2=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,dxmoy,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\langle \Delta x \rangle \ \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig2,'figures/conv0_dx_dx', '-depsc');
    
    fig3=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,pmoy,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\langle p \rangle \ \rm [kg \ m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig3,'figures/conv0_dx_p', '-depsc');

    fig4=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,dpmoy,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\langle \Delta p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig4,'figures/conv0_dx_dp', '-depsc');

    fig5=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2, errx,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0}|} \quad \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig5,'figures/conv0_dx_x_th', '-depsc');

    fig6=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2, errp,'k+');
    xlabel('$(\Delta x)^2 \ \rm [m^2]$','Interpreter','Latex')
    ylabel('$\max_t{|\frac{\langle p \rangle - p_{\rm th}}{A_0 \omega}|} \quad \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('Valeurs numériques', 'Régression linéaire');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig6,'figures/conv0_dx_p_th', '-depsc');
    
    
elseif(strcmp(paraName,'delta') || strcmp(paraName,'V0'))
    
    fig1=figure('Position',[50,50,600,450]);
    hold on
    h=plot(V0, prob_dmax,[E0 E0],[0 1],'r--');
    hold off
    xlabel('$V_0 \ \rm [J]$','Interpreter','Latex')
    ylabel('$\max_t{P_{x>0}}$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
%     lgd=legend('show');
%     set(lgd,'fontsize',14,'Location','northwest');
    print(fig1,'figures/delta', '-depsc');
    
    if nsimul<6
        fig2=figure('Position',[50,50,600,450]);
        hold on
        for i=1:nsimul
            p=plot(t, prob_g{i} ,'--','DisplayName',sprintf('$V_0 = %.1fE_0$',V0(i)/E0));
            h(i)=plot(t, prob_d{i} ,'-','DisplayName',sprintf('$V_0 = %.1fE_0$',V0(i)/E0),'Color',get(p,'Color'));
        end
        hold off
        xlabel('$t \ \rm [s]$','Interpreter','Latex')
        ylabel("Probabilit\'{e}",'Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
        set(gca,'FontSize',25)
        grid on
        lgd=legend(h,'Interpreter','Latex');
        set(lgd,'fontsize',14,'Location','southeast');
        print(fig2,'figures/delta', '-depsc');
        
        for i=1:nsimul
            fig(i)=figure('Position',[50,50,650,450]);
            pcolor(x,t,psi2{i})
            shading interp
            colormap jet
            c = colorbar;
            xlabel('$x \ \rm [m]$','Interpreter','Latex')
            ylabel('$t \ \rm [s]$','Interpreter','Latex')
            ylabel(c,'$|\psi(x,t)|^2$','Interpreter','Latex')
            grid on, box on
            set(gca,'FontSize',22)
            title(sprintf('$V_0 = %.1fE_0$',V0(i)/E0),'Interpreter','Latex')
            pos=get(gca,'position');  % retrieve the current values
            pos(3)=0.9*pos(3);        % try reducing width 10%
            set(gca,'position',pos);  % write the new values
            print(fig(i),sprintf('figures/V0_%d',i), '-depsc');
        end
    end
    
elseif(strcmp(paraName,'n'))
    
    fig1=figure('Position',[50,50,600,450]);
    hold on
    l=round(3*length(t)/5);
    for i=1:nsimul
        h=plot(t(1:l), prob_d{i}(1:l),'DisplayName',sprintf('n = %d',n(i)));
    end
    hold off
    xlabel('$t \ \rm [s] $','Interpreter','Latex')
    ylabel('$P_{x>0}$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
    set(gca,'FontSize',25)
    set(h,'MarkerSize',11)
    grid on
    ylim([0 1])
    lgd=legend('show');
    set(lgd,'fontsize',14,'Location','northeast');
    print(fig1,'figures/n_scan', '-depsc');
   
end