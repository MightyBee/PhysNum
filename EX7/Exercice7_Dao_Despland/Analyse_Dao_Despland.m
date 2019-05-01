%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fichier = 'output';

% excitation='harmo';
excitation='pulse';

% bord='fixe';
% bord='libre';
bord='sortie';

% fichier = [excitation '_' bord];
fichier = 'cfl1';

data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);
[X,T]=meshgrid(x,t);

%%
h0=u.^2/9.81;

%% Figures %%
%%%%%%%%%%%%%
figure
plot(x,f(end,:));
hold on 
plot(x,9*h0.^(-1/4),'--')
grid
xlabel('x [m]')
ylabel('f(x,t_{fin}) [m]')


fig1=figure('Position',[50,50,600,450]);
pcolor(x,t,f)
shading interp
colormap jet
c = colorbar;
xlabel('$x \ \rm [m]$','Interpreter','Latex')
ylabel('$t \ \rm [s]$','Interpreter','Latex')
ylabel(c,'$f(x,t) \ \rm [m]$','Interpreter','Latex')

title('$\beta_{\rm CFL} = 0.1$','Interpreter','Latex')
grid on, box on
set(gca,'FontSize',25)
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
print(fig1,sprintf('figures/reflexion_%s',fichier), '-depsc');
% print(fig1,sprintf('figures/%s',fichier), '-depsc');

fig2=figure('Position',[50,50,600,450]);
semilogy(t,E)
xlabel('$t \ \rm [s]$','Interpreter','Latex')
ylabel('$E \ \rm [J]$','Interpreter','Latex')
grid on
set(gca,'FontSize',25)
title('$\beta_{\rm CFL} = 0.1$','Interpreter','Latex')%,'FontSize',8)
% pos=get(gca,'position');  % retrieve the current values
% pos(3)=0.9*pos(3);        % try reducing width 10%
% set(gca,'position',pos);  % write the new values
print(fig2,sprintf('figures/energie_%s',fichier), '-depsc');
% print(fig1,sprintf('figures/%s',fichier), '-depsc');

figure('Name',['Analyse de ' fichier])
subplot(2,2,1)
plot(x,u)
grid
xlabel('x [m]')
ylabel('u [m/s]')

subplot(2,2,2)
plot(t,E)
grid
xlabel('t [s]')
ylabel('E [m^3]')

subplot(2,2,4)
pcolor(x,t,f)
shading interp
colormap jet
c = colorbar;
xlabel('x [m]')
ylabel('t [s]')
ylabel(c,'f(x,t) [m]')
% s=surf(X,T,f);
% grid
% xlabel('x [m]')
% ylabel('t [s]')
% zlabel('f(x,t) [m]')
% zlim([min(f(:)),max(f(:))])
% s.EdgeColor = 'none';

subplot(2,2,3)
h = plot(x,f(1,:));
grid
xlabel('x [m]')
ylabel('f(x,t) [m]')
ht = title('t=0 s');
ylim([min(f(:)),max(f(:))])

% w=waitforbuttonpress;
for i=2:length(t)
    pause(.01)
    if ~ishandle(h)
         break % Arrete l'animation si la fenetre est fermee
     end
     set(h,'YData',f(i,:))
     set(ht,'String',sprintf('t=%0.2f s',t(i)))
end
 
