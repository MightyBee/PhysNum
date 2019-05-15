%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier='detection_N';

data = load([fichier,'_obs.out']);
t_N      = data(:,1);
prob_g_N = data(:,2);
prob_d_N = data(:,3);

fichier='detection_Y';
data = load([fichier,'_obs.out']);
t      = data(:,1);
prob_g = data(:,2);
prob_d = data(:,3);
E      = data(:,4);
xmoy   = data(:,5);
x2moy  = data(:,6);
pmoy   = data(:,7);
p2moy  = data(:,8);

data = load([fichier,'_pot.out']);
x = data(:,1);
V = data(:,2);

psi2 = load([fichier,'_psi2.out']);
psi = sqrt(psi2);

omega=0.003;


%% Analyse %%
%%%%%%%%%%%%%
dxmoy=sqrt(x2moy-xmoy.^2);
dpmoy=sqrt(p2moy-pmoy.^2);

%% Figures %%
%%%%%%%%%%%%%
%%{
%%
fig1=figure('Position',[50,50,600,450]);
hold on
p=plot(t, prob_g_N ,'--');
h(1)=plot(t, prob_d_N ,'-','DisplayName',"Sans d\'{e}tection",'Color',get(p,'Color'));
p=plot(t, prob_g ,'--','DisplayName',sprintf('$V_0 = %.1fE_0$',V0(i)/E0));
h(2)=plot(t, prob_d ,'-','DisplayName',"Avec d\'{e}tection",'Color',get(p,'Color'));
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel("Probabilit\'{e}",'Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
set(gca,'FontSize',25)
grid on
lgd=legend(h,'Interpreter','Latex');
set(lgd,'fontsize',17,'Location','northeast');
print(fig1,'figures/prob_detect', '-depsc');

%%
fig3=figure('Position',[50,50,600,450]);
hold on
plot(t,E)
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$E \ \rm [J]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
ylim([0 1.1*max(abs(E))])
set(gca,'FontSize',25)
grid on
print(fig3,'figures/E_detect', '-depsc');

%%
fig4=figure('Position',[50,50,600,450]);
xlim([0 5000])
hold on
plot(t,xmoy)
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle x \rangle \ \rm [m]$','Interpreter','Latex')
set(gca,'FontSize',22)
grid on
print(fig4,'figures/x_detect', '-depsc');

%%
fig5=figure('Position',[50,50,600,450]);
plot(t,pmoy,'DisplayName',"Th\'eorie quantique")
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex')
xlim([0 5000])
set(gca,'FontSize',22)
grid on
print(fig5,'figures/p_detect', '-depsc');

%%
fig6=figure('Position',[50,50,600,450]);
plot(t,dxmoy)
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle \Delta x \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') 
xlim([0 5000])
set(gca,'FontSize',22)
grid on
print(fig6,'figures/dx_detect', '-depsc');

%%
fig7=figure('Position',[50,50,600,450]);
plot(t,dpmoy)
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle \Delta p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex')
xlim([0 5000])
set(gca,'FontSize',22)
grid on
print(fig7,'figures/dp_detect', '-depsc');

%%
fig8=figure('Position',[50,50,650,450]);
pcolor(x,t,psi2)
shading interp
colormap jet
c = colorbar;
xlabel('$x \ \rm [m]$','Interpreter','Latex')
ylabel('$t \ \rm [s]$','Interpreter','Latex')
ylabel(c,'$|\psi(x,t)|^2$','Interpreter','Latex')
grid on, box on
set(gca,'FontSize',22)
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
print(fig8,'figures/psi2__detect', '-depsc');

%}

%%{
%%
figure
yyaxis left
h = plot(x,psi(1,:));
ylabel('\psi(x,t) [m]')
ylim([min(psi(:)),max(psi(:))])
yyaxis right
plot(x,V,'--')
grid
xlabel('x [m]')
ylabel('V(x) [J]')
ht = title('t=0 s');


w=waitforbuttonpress;
stop=1;

for i=2:1:length(t)
    pause(.01)
    if ~ishandle(h)
        break % Arrete l'animation si la fenetre est fermee
    end
    if t(i)>=1000 && stop
        w=waitforbuttonpress;
        stop=0;
        set(h,'YData',psi(i,:))
        set(ht,'String',sprintf('t=%0.2f s',t(i)))
        w=waitforbuttonpress;
    end
    set(h,'YData',psi(i,:))
    set(ht,'String',sprintf('t=%0.2f s',t(i)))
end
%}
