%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output_64';

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


P0=2*pi*14/400;
E0=P0^2/2; %mean(data(:,4));
A=sqrt(2*E0)/omega;
A=max(xmoy);
fprintf('\n E=%.15g \n\n',E(1))


%% Figures %%
%%%%%%%%%%%%%
%%{
%%
fig1=figure('Position',[50,50,600,450]);
hold on
plot(t,prob_g,'DisplayName','$ P_{x<0} $')
plot(t,prob_d,'DisplayName','$ P_{x>0} $')
plot(t,prob_g+prob_d,'--','DisplayName','$ P_{x<0}+P_{x>0} $')
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel("Probabilit\'{e}",'Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
ylim([0 1.1])
set(gca,'FontSize',25)
grid on
lgd=legend('show','Interpreter','Latex');
set(lgd,'fontsize',17,'Location','southeast');
print(fig1,'figures/prob', '-depsc');

%%
fig1b=figure('Position',[50,50,600,450]);
plot(t,abs(prob_g+prob_d-1),'-','DisplayName','$ P_{x<0}+P_{x>0} $')
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel("$|P_{x<0}+P_{x>0}-1|$",'Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
% ylim([0 1.1])
set(gca,'FontSize',23)
ax = gca;
ax.YAxis.Exponent = -14;
grid on
% lgd=legend('show','Interpreter','Latex');
% set(lgd,'fontsize',17,'Location','southeast');
print(fig1b,'figures/prob_zoom', '-depsc');

%%
fig2=figure('Position',[50,50,600,450]);
hold on
plot(t,E-E(1))
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$E(t)-E(0) \ \rm [J]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
% ylim([0 1.1])
set(gca,'FontSize',25)
grid on
print(fig2,'figures/E_diff', '-depsc');

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
print(fig3,'figures/E_abs', '-depsc');

%%
fig4=figure('Position',[50,50,600,450]);
hold on
plot(t,dxmoy.*dpmoy,'DisplayName','$ \langle \Delta x \rangle \cdot \langle \Delta p \rangle $')
plot([t(1) t(end)],[0.5 0.5],'--','DisplayName','$ \hbar/2 $')
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle \Delta x \rangle \cdot \langle \Delta p \rangle  \ \rm [kg \, m^2/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
yl=ylim;
ylim([0.4 yl(2)])
set(gca,'FontSize',22)
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
grid on
lgd=legend('show','Interpreter','Latex');
set(lgd,'fontsize',17,'Location','southeast');
print(fig4,'figures/dx_dp', '-depsc');

%%
fig4=figure('Position',[50,50,600,450]);
hold on
plot(t,xmoy,'DisplayName',"Th\'eorie quantique")
plot(t,A*sin(omega*t),'--','DisplayName',"Th\'eorie classique")
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle x \rangle \ \rm [m]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
% ylim([0.4 0.8])
set(gca,'FontSize',22)
% pos=get(gca,'position');  % retrieve the current values
% pos(3)=0.9*pos(3);        % try reducing width 10%
% set(gca,'position',pos);  % write the new values
grid on
lgd=legend('show','Interpreter','Latex');
set(lgd,'fontsize',17,'Location','southeast');
print(fig4,'figures/x_th', '-depsc');

%%
fig5=figure('Position',[50,50,600,450]);
hold on
plot(t,pmoy,'DisplayName',"Th\'eorie quantique")
plot(t,A*omega*cos(omega*t),'--','DisplayName',"Th\'eorie classique")
hold off
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
% ylim([0.4 0.8])
set(gca,'FontSize',22)
% pos=get(gca,'position');  % retrieve the current values
% pos(3)=0.9*pos(3);        % try reducing width 10%
% set(gca,'position',pos);  % write the new values
grid on
lgd=legend('show','Interpreter','Latex');
set(lgd,'fontsize',17,'Location','southeast');
print(fig5,'figures/p_th', '-depsc');

%%
fig6=figure('Position',[50,50,600,450]);
plot(t,dxmoy)
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle \Delta x \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
set(gca,'FontSize',22)
grid on
print(fig6,'figures/dx_moy', '-depsc');

%%
fig7=figure('Position',[50,50,600,450]);
plot(t,dpmoy)
xlabel('$ t \ \rm [s]$','Interpreter','Latex')
ylabel('$\langle \Delta p \rangle \ \rm [kg \, m/s]$','Interpreter','Latex') %$f(5{\rm m},1.5{\rm s}) \ \rm [m]$
xlim([0 5000])
set(gca,'FontSize',22)
grid on
print(fig7,'figures/dp_moy', '-depsc');

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
print(fig8,sprintf('figures/psi2_%s',fichier), '-depsc');

%%
fig9=figure('Position',[50,50,750,400]);
plot(x,V)
set(gca,'FontSize',14)
xlabel('$ x \ \rm [m]$','Interpreter','Latex','FontSize',22)
ylabel('$ V \ \rm [J]$','Interpreter','Latex','FontSize',22)
grid on
print(fig9,'figures/potentiel', '-depsc');

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
