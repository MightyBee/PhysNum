%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = '3';
data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);

%% Figures %%
%%%%%%%%%%%%%
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

subplot(2,2,3)
h = plot(x,f(1,:));
grid
xlabel('x [m]')
ylabel('f(x,t) [m]')
ht = title('t=0 s');
ylim([min(f(:)),max(f(:))])
% for i=2:length(t)
%     pause(.01)
%     if ~ishandle(h)
%         break % Arrete l'animation si la fenetre est fermee
%     end
%     set(h,'YData',f(i,:))
%     set(ht,'String',sprintf('t=%0.2f s',t(i)))
% end
