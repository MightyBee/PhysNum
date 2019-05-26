set(groot,'defaultAxesFontSize',15)
set(groot,'defaultAxesLabelFontSizeMultiplier',22/15)
set(0,'defaultTextInterpreter','latex');
set(0,'defaultLegendInterpreter','latex');

%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier='detect_0';

data = load([fichier,'_obs.out']);
t_N      = data(:,1);
prob_g_N = data(:,2);
prob_d_N = data(:,3);
E_N      = data(:,4);
xmoy_N   = data(:,5);
x2moy_N  = data(:,6);
pmoy_N   = data(:,7);
p2moy_N  = data(:,8);

fichier='detect_1';
data = load([fichier,'_obs.out']);
t_1      = data(:,1);
prob_g_1 = data(:,2);
prob_d_1 = data(:,3);
E_1      = data(:,4);
xmoy_1   = data(:,5);
x2moy_1  = data(:,6);
pmoy_1   = data(:,7);
p2moy_1  = data(:,8);

fichier='detect_3';
data = load([fichier,'_obs.out']);
t_2      = data(:,1);
prob_g_2 = data(:,2);
prob_d_2 = data(:,3);
E_2      = data(:,4);
xmoy_2   = data(:,5);
x2moy_2  = data(:,6);
pmoy_2   = data(:,7);
p2moy_2  = data(:,8);

data = load([fichier,'_pot.out']);
x = data(:,1);
V = data(:,2);

psi2 = load([fichier,'_psi2.out']);
psi = sqrt(psi2);

omega=0.003;


%% Analyse %%
%%%%%%%%%%%%%
dxmoy_N=sqrt(x2moy_N-xmoy_N.^2);
dpmoy_N=sqrt(p2moy_N-pmoy_N.^2);
dxmoy=sqrt(x2moy_1-xmoy_1.^2);
dpmoy=sqrt(p2moy_1-pmoy_1.^2);

%% Figures %%
%%%%%%%%%%%%%
%%{
%%
fig1=figure('Position',[50,50,600,450]);
hold on
% p=plot(t_1, prob_g_N ,'--');
% h(2)=plot(t_1, prob_d_N ,'-','DisplayName',"Sans d\'{e}tection",'Color',get(p,'Color'));
% p=plot(t_1, prob_g_1 ,'--');
% h(1)=plot(t_1, prob_d_1 ,'-','DisplayName',"Avec d\'{e}tection",'Color',get(p,'Color'));
% p=plot(t_2, prob_g_2 ,'--');
% h(1)=plot(t_2, prob_d_2 ,'-','DisplayName',"Avec d\'{e}tection",'Color',get(p,'Color'));
h(1)=plot(t_1, prob_d_N ,'--','DisplayName',"Sans d\'{e}tection");
h(2)=plot(t_1, prob_d_1 ,'-.','DisplayName',"D\'{e}tection \`a droite");
h(3)=plot(t_2, prob_d_2 ,'-','DisplayName',"D\'{e}tection \`a gauche");
hold off
yticks([0 0.25 0.5 0.75 1])
xlabel('$t$')
ylabel("$P_{x>0}$")
xlim([0 5000])
grid on, box on
lgd=legend(h);
set(lgd,'Location','northeast');
print(fig1,'figures/prob_detect', '-depsc');

%%
fig3=figure('Position',[50,50,600,450]);
hold on
h(2)=plot(t_1,E_N,'--','DisplayName',"Sans d\'{e}tection");
h(1)=plot(t_1,E_1,'DisplayName',"Avec d\'{e}tection");
hold off
xlabel('$t$')
ylabel('$E$') 
xlim([0 5000])
ylim([0 1.1*max(abs(E_1))])
grid on, box on
lgd=legend(h);
set(lgd,'Location','southeast');
print(fig3,'figures/E_detect', '-depsc');

%%
fig4=figure('Position',[50,50,600,450]);
hold on
h(1)=plot(t_1,xmoy_N,'--','DisplayName',"Sans d\'{e}tection");
h(2)=plot(t_1,xmoy_1,'-.','DisplayName',"D\'{e}tection \`a droite");
h(3)=plot(t_1,xmoy_2,'-','DisplayName',"D\'{e}tection \`a gauche");
hold off
xlim([0 5000])
xlabel('$t$')
ylabel('$\langle x \rangle$')
grid on, box on
lgd=legend(h);
set(lgd,'Location','northeast');
print(fig4,'figures/x_detect', '-depsc');

%%
fig5=figure('Position',[50,50,600,450]);
hold on
h(2)=plot(t_1,pmoy_N,'--','DisplayName',"Sans d\'{e}tection");
h(1)=plot(t_1,pmoy_1,'DisplayName',"Avec d\'{e}tection");
hold off
xlabel('$t$')
ylabel('$\langle p \rangle$')
xlim([0 5000])
grid on, box on
lgd=legend(h);
set(lgd,'Location','northeast');
print(fig5,'figures/p_detect', '-depsc');

%%
fig6=figure('Position',[50,50,600,450]);
hold on
h(2)=plot(t_1,dxmoy_N,'--','DisplayName',"Sans d\'{e}tection");
h(1)=plot(t_1,dxmoy,'DisplayName',"Avec d\'{e}tection");
hold off
xlabel('$t$')
ylabel('$\langle \Delta x \rangle$') 
xlim([0 5000])
grid on, box on
lgd=legend(h);
set(lgd,'Location','southeast');
print(fig6,'figures/dx_detect', '-depsc');

%%
fig7=figure('Position',[50,50,600,450]);
hold on
h(2)=plot(t_1,dpmoy_N,'--','DisplayName',"Sans d\'{e}tection");
h(1)=plot(t_1,dpmoy,'DisplayName',"Avec d\'{e}tection");
hold off
xlabel('$t$')
ylabel('$\langle \Delta p \rangle$')
xlim([0 5000])
grid on, box on
lgd=legend(h);
set(lgd,'Location','southeast');
print(fig7,'figures/dp_detect', '-depsc');

%%
fig8=figure('Position',[50,50,650,450]);
pcolor(x,t_1,psi2)
shading interp
colormap jet
c = colorbar;
xlabel('$x$')
ylabel('$t$')
ylabel(c,'$|\psi(x,t)|^2$')
grid on, box on
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
print(fig8,'figures/psi2_detect', '-depsc');

%}

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
stop=1;

for i=2:1:length(t_1)
    pause(.01)
    if ~ishandle(h)
        break % Arrete l'animation si la fenetre est fermee
    end
    if t_1(i)>=1000 && stop
        w=waitforbuttonpress;
        stop=0;
        set(h,'YData',psi(i,:))
        set(ht,'String',sprintf('$t=%0.2f$',t_1(i)))
        w=waitforbuttonpress;
    end
    set(h,'YData',psi(i,:))
    set(ht,'String',sprintf('$t=%0.2f$',t_1(i)))
end
%}
