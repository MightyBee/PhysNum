set(groot,'defaultAxesFontSize',15)
set(groot,'defaultAxesLabelFontSizeMultiplier',22/15)
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

fichier = 'quasi_stat';

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
E0=P0^2/2;
E_moy=mean(E);
A=sqrt(2*E_moy)/omega;
% A=max(xmoy);
fprintf('\n E=%.15g \n\n',E(1))


%% Figures %%
%%%%%%%%%%%%%
%%
fig1=figure('Position',[50,50,600,450]);
hold on
plot(t,prob_g,'DisplayName','$ P_{x<0} $')
plot(t,prob_d,'DisplayName','$ P_{x>0} $')
plot(t,prob_g+prob_d,'--','DisplayName','$ P_{x<0}+P_{x>0} $')
hold off
set(gca,'FontSize',15)
xlabel('$ t $','Interpreter','Latex','FontSize',22)
ylabel("Probabilit\'{e}",'Interpreter','Latex','FontSize',22)
xlim([0 5000])
ylim([0 1.1])
grid on
lgd=legend('show','Interpreter','Latex');
set(lgd,'fontsize',15,'Location','southeast');
print(fig1,'figures/prob', '-depsc');

%%
fig3=figure('Position',[50,50,600,450]);
plot(t,E)
set(gca,'FontSize',15)
xlabel('$ t $','Interpreter','Latex','FontSize',22)
ylabel('$E $','Interpreter','Latex','FontSize',22)
xlim([0 5000])
ylim([0 1.1*max(abs(E))])
grid on
print(fig3,'figures/E_abs', '-depsc');

%%
fig4=figure('Position',[50,50,600,450]);
plot(t,xmoy)
set(gca,'FontSize',15)
xlabel('Temps','Interpreter','Latex','FontSize',22)
ylabel('Position','Interpreter','Latex','FontSize',22)
xlim([0 5000])
grid on
print(fig4,'figures/x_th', '-depsc');

%%
fig5=figure('Position',[50,50,600,450]);
hold on
plot(t,pmoy)
set(gca,'FontSize',15)
xlabel('Temps','Interpreter','Latex','FontSize',22)
ylabel("Quantit\'e de mouvement",'Interpreter','Latex','FontSize',22)
xlim([0 5000])
grid on
print(fig5,'figures/p_th', '-depsc');

%%
fig6=figure('Position',[50,50,600,450]);
plot(t,dxmoy)
set(gca,'FontSize',15)
xlabel('$ t $','Interpreter','Latex','FontSize',22)
ylabel('$\langle \Delta x \rangle $','Interpreter','Latex','FontSize',22)
xlim([0 5000])
grid on
print(fig6,'figures/dx_moy', '-depsc');

%%
fig8=figure('Position',[50,50,650,450]);
pcolor(x,t,psi2)
set(gca,'FontSize',16)
shading interp
colormap jet
c = colorbar;
xlabel('$x$','Interpreter','Latex','FontSize',22)
ylabel('$t$','Interpreter','Latex','FontSize',22)
ylabel(c,'$|\psi(x,t)|^2$','Interpreter','Latex','FontSize',22)
grid on, box on
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
print(fig8,sprintf('figures/psi2_%s',fichier), '-depsc');

%%
fig9=figure('Position',[50,50,650/59*40,400]);
plot(x,V)
set(gca,'FontSize',16)
xlabel('$ x $','Interpreter','Latex','FontSize',22)
ylabel('$ V $','Interpreter','Latex','FontSize',22)
grid on
print(fig9,sprintf('figures/potentiel_%s',fichier), '-depsc');

%%
fig10=figure('Position',[50,50,700,700]);

subplot('Position',[0.16,0.1,0.6,0.5])
pcolor(x,t,psi2)
set(gca,'FontSize',14)
shading interp
colormap jet
c = colorbar;
xlabel('$x$','Interpreter','Latex','FontSize',22)
ylabel('$t$','Interpreter','Latex','FontSize',22)
ylabel(c,'$|\psi(x,t)|^2$','Interpreter','Latex','FontSize',22)
grid on, box on
pos=get(gca,'position');  % retrieve the current values
pos(3)=1.225*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
subplot('Position',[0.16,0.69,0.6,0.28])
plot(x,V)
set(gca,'FontSize',14)
xlabel('$ x $','Interpreter','Latex','FontSize',22)
ylabel('$ V $','Interpreter','Latex','FontSize',22)
grid on, box off
print(fig10,sprintf('figures/duo_%s',fichier), '-depsc');

%%{
%%
figure
yyaxis left
h = plot(x,psi(1,:));
ylabel('$\psi(x,t)$')
ylim([min(psi(:)),max(psi(:))])
yyaxis right
plot(x,V,'--')
grid
xlabel('$x$')
ylabel('$V(x)$')
ht = title('$t=0 $');

w=waitforbuttonpress;

for i=2:3:length(t)
    pause(.01)
    if ~ishandle(h)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(h,'YData',psi(i,:))
    set(ht,'String',sprintf('$t=%0.2f$',t(i)))
end
%}

% ca=input('Tout fermer [y/n]  ? ','s');
% if strcmp(ca,'y')
%     close all;
% end
