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

nsimul = 151; % Nombre de simulations a faire
paraName='delta';

m=1;
omega=0.003;
P0=2*pi*14/400;
E0=P0^2/2; %mean(data(:,4));


if strcmp(paraName,'dt')
    N=round(logspace(3,4, nsimul));
    dt=5000./N;
    paramstr = {"dt"}; % Nom du parametre a scanner
    param = dt; % Valeurs du parametre a scanner
elseif strcmp(paraName,'Ninters')
    Ninters=round(logspace(1.9,3.1,nsimul));
    paramstr = {"Ninters"}; % Nom du parametre a scanner
    param = [Ninters]; % Valeurs du parametre a scanner
elseif strcmp(paraName,'delta')
%     delta=linspace(0,150,nsimul);
    V0=linspace(0,2.5,nsimul)*E0;
    delta=sqrt(2*V0/(m*omega^2));
    paramstr = {"delta"; "x0"}; % Nom du parametre a scanner
    param = [delta; -delta]; % Valeurs du parametre a scanner
elseif strcmp(paraName,'V0')
    V0=[0.2 1 2]*E0;
    nsimul=length(V0);
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
elseif strcmp(paraName,'n_cont')
    n=linspace(11,12,nsimul);
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
%     system(cmd);
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
    E_m = zeros(1,nsimul);
elseif(strcmp(paraName,'V0'))
    prob_g = cell(1,nsimul);
    prob_d = cell(1,nsimul);
    prob_dmax = zeros(1,nsimul);
    psi2 = cell(1,nsimul);
    V=cell(1,nsimul);
elseif(strcmp(paraName,'n') || strcmp(paraName,'n_cont'))
    prob_d = cell(1,nsimul);
    prob_dmax = zeros(1,nsimul);
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
        prob_dmax(i)=max(prob_d{i}(1:round(2*length(t)/5)));
        E_m=mean(data(:,4));
    elseif(strcmp(paraName,'V0'))
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_g{i}=data(:,2);
        prob_d{i}=data(:,3);
        prob_dmax(i)=max(data(:,3));
        data = load([output{i},'_pot.out']);
        x = data(:,1);
        V{i}=data(:,2);
        psi2{i} = load([output{i},'_psi2.out']);
    elseif(strcmp(paraName,'n') || strcmp(paraName,'n_cont'))
        data = load([output{i},'_obs.out']);
        t=data(:,1);
        prob_d{i}=data(:,3);
        prob_dmax(i)=max(prob_d{i}(1:round(2*length(t)/5)));
    end
end


%% Figures %%
%%%%%%%%%%%%%

if(strcmp(paraName,'dt'))
    fig1=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,xmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle x \rangle$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig1,'figures/conv0_dt_x', '-depsc');

    fig2=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,dxmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle \Delta x \rangle$','Interpreter','Latex','FontSize',22) 
    grid on
    print(fig2,'figures/conv0_dt_dx', '-depsc');
    
    fig3=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,pmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    lgd=legend('show','Interpreter','Latex');
    set(lgd,'FontSize',15,'Location','northeast');
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle p \rangle$','Interpreter','Latex','FontSize',22) 
    grid on
    print(fig3,'figures/conv0_dt_p', '-depsc');

    fig4=figure('Position',[50,50,600,450]);
    h=plot(dt.^2,dpmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle \Delta p \rangle$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig4,'figures/conv0_dt_dp', '-depsc');

    fig5=figure('Position',[50,50,600,450]);
    h=plot(dt.^2, errx,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0}|}$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig5,'figures/conv0_dt_x_th', '-depsc');

    fig6=figure('Position',[50,50,600,450]);
    h=plot(dt.^2, errp,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta t)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\max_t{|\frac{\langle p \rangle - p_{\rm th}}{A_0 \omega}|}$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig6,'figures/conv0_dt_p_th', '-depsc');
    
   
elseif(strcmp(paraName,'Ninters'))
    
    fig1=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,xmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta h)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle x \rangle$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig1,'figures/conv0_dh_x', '-depsc');

    fig2=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,dxmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta h)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle \Delta x \rangle$','Interpreter','Latex','FontSize',22) 
    grid on
    print(fig2,'figures/conv0_dh_dx', '-depsc');
    
    fig3=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,pmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta h)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle p \rangle$','Interpreter','Latex','FontSize',22) 
    grid on
    print(fig3,'figures/conv0_dh_p', '-depsc');

    fig4=figure('Position',[50,50,600,450]);
    h=plot((400./Ninters).^2,dpmoy,'k+');
    set(gca,'FontSize',15)
    set(h,'MarkerSize',11)
    xlabel('$(\Delta h)^2$','Interpreter','Latex','FontSize',22)
    ylabel('$\langle \Delta p \rangle$','Interpreter','Latex','FontSize',22)
    grid on
    print(fig4,'figures/conv0_dh_dp', '-depsc');

    fig5=figure('Position',[50,50,600,450]);
    h=loglog(400./Ninters, errx,'k+');
    set(gca,'FontSize',18)
    set(h,'MarkerSize',12)
    xlabel('$\Delta h$','Interpreter','Latex','FontSize',26)
    ylabel('$\max_t{|\frac{\langle x \rangle - x_{\rm th}}{A_0}|}$','Interpreter','Latex','FontSize',26)
    xticks([0.5 1 2 3 4 5])
    grid on
    print(fig5,'figures/conv0_dh_x_th', '-depsc');

    fig6=figure('Position',[50,50,600,450]);
    h=loglog(400./Ninters, errp,'k+');
    set(gca,'FontSize',18)
    set(h,'MarkerSize',12)
    xlabel('$\Delta h$','Interpreter','Latex','FontSize',26)
    ylabel('$\max_t{|\frac{\langle p \rangle - p_{\rm th}}{A_0 \omega}|}$','Interpreter','Latex','FontSize',26)
    xticks([0.5 1 2 3 4 5])
    grid on
    print(fig6,'figures/conv0_dh_p_th', '-depsc');

    
    
elseif(strcmp(paraName,'delta') || strcmp(paraName,'V0'))
    
    fig1=figure('Position',[50,50,600,450]);
    hold on
    plot(V0./E_m, prob_dmax,'DisplayName','Quantique');
    plot([0 1 1 2.5],[1 1 0 0],'r--','DisplayName','Classique');
    hold off
    xlim([0 2.5])
    set(gca,'FontSize',15)
    xlabel('$V_0/E_0$','Interpreter','Latex','FontSize',22)
    ylabel('$P_{\rm trans}$','Interpreter','Latex','FontSize',22)
    grid on
    lgd=legend('show','Interpreter','Latex');
    set(lgd,'FontSize',15,'Location','northeast');
    print(fig1,'figures/V0_scan', '-depsc');
    
    if nsimul<6
        fig2=figure('Position',[50,50,600,450]);
        hold on
        for i=1:nsimul
            p=plot(t, prob_g{i} ,'--','DisplayName',sprintf('$V_0 = %.1fE_0$',V0(i)/E0));
            h(i)=plot(t, prob_d{i} ,'-','DisplayName',sprintf('$V_0 = %.1fE_0$',V0(i)/E0),'Color',get(p,'Color'));
        end
        hold off
        set(gca,'FontSize',15)
        xlabel('$t$','Interpreter','Latex','FontSize',22)
        ylabel("Probabilit\'{e}",'Interpreter','Latex','FontSize',22)
        grid on
        lgd=legend(h,'Interpreter','Latex');
        set(lgd,'fontsize',15,'Location','southeast');
        print(fig2,'figures/prob_tunnel', '-depsc');
        
        for i=1:nsimul
            fig(i)=figure('Position',[50,50,650,450]);
            pcolor(x,t,psi2{i})
            shading interp
            colormap jet
            c = colorbar;
            set(gca,'FontSize',15)
            xlabel('$x$','Interpreter','Latex','FontSize',22)
            ylabel('$t$','Interpreter','Latex','FontSize',22)
            ylabel(c,'$|\psi(x,t)|^2$','Interpreter','Latex','FontSize',22)
            grid on, box on
            title(sprintf('$V_0 = %.1fE_0$',V0(i)/E0),'Interpreter','Latex')
            pos=get(gca,'position');  % retrieve the current values
            pos(3)=0.9*pos(3);        % try reducing width 10%
            set(gca,'position',pos);  % write the new values
            print(fig(i),sprintf('figures/V0_%d',i), '-depsc');
            
            fig2(i)=figure('Position',[50,50,600,450]);
            plot(x,V{i})
            set(gca,'FontSize',15)
            xlabel('$x$','Interpreter','Latex','FontSize',22)
            ylabel('$V$','Interpreter','Latex','FontSize',22)
            grid on
            print(fig2(i),sprintf('figures/potentiel_scan%d',i), '-depsc');
        end
    end
    
elseif(strcmp(paraName,'n'))
    
    fig1=figure('Position',[50,50,600,450]);
    hold on
    l=round(3*length(t)/5);
    for i=1:nsimul
        h=plot(t(1:l), prob_d{i}(1:l),'DisplayName',sprintf('$n=%d$',n(i)));
    end
    hold off
    set(gca,'FontSize',15)
    xlabel('$t$','Interpreter','Latex','FontSize',22)
    ylabel('$P_{x>0}$','Interpreter','Latex','FontSize',22)
    set(h,'MarkerSize',11)
    grid on
    ylim([0 1])
    yticks([0 0.25 0.5 0.75 1])
    lgd=legend('show','Interpreter','Latex');
    set(lgd,'fontsize',15,'Location','northeast');
    print(fig1,'figures/n_scan', '-depsc');

elseif(strcmp(paraName,'n_cont'))
   
    fig1=figure('Position',[50,50,600,450]);
    hold on
    plot(n, prob_dmax,'+', 'MarkerSize',11);
    plot(n,0.5*ones(size(n)),'--')
    hold off
    set(gca,'FontSize',16)
    xlabel('$n$','Interpreter','Latex','FontSize',25)
    ylabel('$P_{\rm trans}$','Interpreter','Latex','FontSize',25)
    grid on, box on
    print(fig1,'figures/n_cont_scan', '-depsc');
    
end