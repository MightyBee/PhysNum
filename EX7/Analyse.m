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

    switch(cb_gauche)
    {
      case fixe:
        fnext[0] = fnow[0]; // TODO : Completer la condition au bord gauche fixe
        break;

      case libre:
        fnext[0] = fnext[1]; // TODO : Completer la condition au bord gauche libre
        break;

      case harmonique:
        fnext[0] = A*sin(omega*t); // TODO : Completer la condition au bord gauche harmonique
        break;

      case sortie:
        fnext[0] = fnow[0]+sqrt((*u2)(0.5*dx))*dt/dx*(fnow[1]-fnow[0]); // TODO : Completer la condition au bord gauche "sortie de l'onde"
        break;
    }

    switch(cb_droit)
    {
      case fixe:
      fnext[N-1] = fnow[N-1]; // TODO : Completer la condition au bord droit fixe
        break;

      case libre:
        fnext[N-1] = fnext[N-2]; // TODO : Completer la condition au bord droit libre
        break;

      case harmonique:
        fnext[N-1] = A*sin(omega*t); // TODO : Completer la condition au bord droit harmonique
        break;

      case sortie:
        fnext[N-1] = fnow[0]-sqrt((*u2)((N-0.5)*dx))*dt/dx*(fnow[N-1]-fnow[N-2]);; // TODO : Completer la condition au bord droit "sortie de l'onde"
        break;
    }
    
    
    double E(0.0);
  for(size_t i(0); i<f.size()-1; i++){
    E+=(f[i]+f[i+1])*0.5*dx;
  }
  return E;
