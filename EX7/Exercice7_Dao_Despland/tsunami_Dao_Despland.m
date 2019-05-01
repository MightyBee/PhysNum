%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration_t.in'; % Nom du fichier d'entree de base
dossier='simulations';


N=2330;
liste='ABC';

x=cell(3,1);
u=cell(3,1);
x_=cell(3,1);
v_=cell(3,1);
h_th=cell(3,1);
peak_=cell(3,1);

%% Simulations %%
%%%%%%%%%%%%%%%%%

for num=1:3
    schema=liste(num);
    output = sprintf('%s/schema=%s_N=%d',dossier,schema,N);
    cmd = sprintf('%s%s %s cb_gauche=harmonique Npoints=%d schema=%s output=%s', repertoire, executable, input, N, schema, output);
    disp(cmd)
    % system(cmd);
    
    %% Chargement des resultats %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fichier = output;
    data = load([fichier,'_u.out']);
    x{num} = data(:,1);
    u{num} = data(:,2);
    data = load([fichier,'_E.out']);
    t = data(:,1);
    E = data(:,2);
    data = load([fichier,'_f.out']);
    f = data(:,2:end);
    
    %% Analyse %%
    %%%%%%%%%%%%%
    m=100; % nb points ou la vitesse est evaluee
    n=40; % nb de pas spacial pour l'evaluation de la vitesse
    dx=max(x{num}(2:end)-x{num}(1:end-1));
    h0=u{num}.^2/9.81;
    
    v_{num}=zeros(m,1);
    u_=zeros(m,1);
    x_{num}=zeros(m,1);
    peak_{num}=zeros(m,1);
    j=1;
    
    for a=round(linspace(2*n,length(x{num})-2*n,m))
        
        [peak1,t1]=findpeaks(f(:,a-n),t,'MinPeakDistance',800);
        i=find(t==t1(1));
        [~,t1_]=inter_max(t(i-4:i+4),f(i-4:i+4,a-n),2);
        
        [peak,t1]=findpeaks(f(:,a),t,'MinPeakDistance',800);
        peak_{num}(j)=peak(3);
        
        [peak2,t2]=findpeaks(f(:,a+n),t,'MinPeakDistance',800);
        i=find(t==t2(1));
        [~,t2_]=inter_max(t(i-4:i+4),f(i-4:i+4,a+n),3);
        
        v_{num}(j)=(x{num}(a+n)-x{num}(a-n))/(t2_-t1_);
        u_(j)=u{num}(a);
        x_{num}(j)=x{num}(a);
        
        j=j+1;
        
    end
    
    err=mean(abs(v_{num}-u_));
    i1=round(400/800*length(h0));
    i2=round(400/800*m);
    disp(h0.^(1/4));
    disp(peak_{num}(i2));
    disp(h0(i1)^(1/4));
    if strcmp(schema,'A')
        h_th{num}=(h0./h0(i1)).^(1/4)*peak_{num}(i2);
    elseif strcmp(schema,'B')
        h_th{num}=h0.^(-1/4)*peak_{num}(i2)/h0(i1)^(-1/4);
        %    h_th{num}=h_th{num}*h_th{num}(i0)^(-1/4)/h0(i0)^(-1/4);
    elseif strcmp(schema,'C')
        h_th{num}=h0.^(-3/4)*peak_{num}(i2)/h0(i1)^(-3/4);
        %    h_th{num}=h_th{num}*h_th{num}(i0)^(-3/4)/h0(i0)^(-3/4);
    end
end
%% Figures %%
%%%%%%%%%%%%%

% end
% figure
% plot(t,f(:,a),'LineWidth',1)
% hold on;
% plot(t1,peak1,'o')
% hold off
% xlabel('Period $T$ [s]')
% ylabel('Amplitude $A$ [V]')
% grid on, box on
%
fig1=figure('Position',[50,50,600,450]);
plot(x{1},u{1})
hold on
plot(x_{1},v_{1},'+',x_{2},v_{2},'+',x_{3},v_{3},'+')
hold off
xlabel('$x \ \rm [m]$','Interpreter','Latex')
ylabel('$u \ \rm [m/s]$','Interpreter','Latex')
grid on, box on
set(gca,'FontSize',25)
%     set(h,'MarkerSize',12)
grid on
lgd=legend('Théorique', 'Schéma A', 'Schéma B', 'Schéma C');
set(lgd,'fontsize',14,'Location','southeast');
print(fig1,sprintf('figures/vitesse_%s',schema), '-depsc');



fig2=figure('Position',[50,50,600,450]);
hold on
h1=plot(x{1},h_th{1});
plot(x_{1},peak_{1},'+','Color',get(h1,'Color'))
h2=plot(x{2},h_th{2});
plot(x_{2},peak_{2},'+','Color',get(h2,'Color'))
h3=plot(x{3},h_th{3});
plot(x_{3},peak_{3},'+','Color',get(h3,'Color'))
hold off
xlabel('$x \ \rm [m]$','Interpreter','Latex')
ylabel('$\hat{f} \ \rm [m]$','Interpreter','Latex')
grid on, box on
set(gca,'FontSize',25)
%     set(h,'MarkerSize',12)
grid on
% ylim([0 5])
lgd=legend('$\propto h(x)^{1/4}$','Schéma A','$\propto h(x)^{-1/4}$', 'Schéma B','$\propto h(x)^{-3/4}$', 'Schéma C');
set(lgd,'fontsize',14,'Location','northeast','Interpreter','Latex');
print(fig2,sprintf('figures/hauteur'), '-depsc');




