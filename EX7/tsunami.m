%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration_t.in'; % Nom du fichier d'entree de base
dossier='simulations/';

nsimul = 10; % Nombre de simulations a faire


% N = round((logspace(1,3, nsimul)./4))*4+1;
N=round(logspace(3,3.7, nsimul));

paramstr = 'Npoints'; % Nom du parametre a scanner
param = N; % Valeurs du parametre a scanner


%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
% for i = 1:nsimul
i=7;
    output{i} = [dossier,paramstr, '=', num2str(param(i))];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%d schema=B output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
% end


%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
err=zeros(nsimul,1);
% for gg=1:nsimul
    gg=7;
fichier = output{gg};
data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_E.out']);
t = data(:,1);
E = data(:,2);
data = load([fichier,'_f.out']);
f = data(:,2:end);

%% Analyse %%
%%%%%%%%%%%%%
m=100;
x0=600000;
n=40;
dx=max(x(2:end)-x(1:end-1));
a = floor(x0/dx+1);
i_ = round(linspace(2*n,length(x)-2*n,m));
v_=zeros(m,1);
u_=zeros(m,1);
x_=zeros(m,1);
j=1;
for a=i_

[peak1,t1]=findpeaks(f(:,a-n),t,'MinPeakDistance',800);
i=find(t==t1(1));
[v,t1_]=inter_max(t(i-4:i+4),f(i-4:i+4,a-n),2);

[peak2,t2]=findpeaks(f(:,a+n),t,'MinPeakDistance',800);
i=find(t==t2(1));
[v,t2_]=inter_max(t(i-4:i+4),f(i-4:i+4,a+n),3);

v_(j)=(x(a+n)-x(a-n))/(t2_-t1_);
u_(j)=u(a);
x_(j)=x(a);
j=j+1;
end

err(gg)=mean(abs(v_-u_))

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
figure
plot(x,u)
hold on
plot(x_,v_,'+')
hold off
xlabel('Period $T$ [s]')
ylabel('Amplitude $A$ [V]')
grid on, box on

figure
plot(N,err)
hold off
xlabel('Period $T$ [s]')
ylabel('Amplitude $A$ [V]')
grid on, box on


% figure('Name',['Analyse de ' fichier])
% subplot(2,2,1)
% plot(x,u)
% grid
% xlabel('x [m]')
% ylabel('u [m/s]')
% 
% subplot(2,2,2)
% plot(t,E)
% grid
% xlabel('t [s]')
% ylabel('$E \rm [m^3]$')
% 
% subplot(2,2,4)
% pcolor(x,t,f)
% shading interp
% colormap jet
% c = colorbar;
% xlabel('x [m]')
% ylabel('t [s]')
% ylabel(c,'f(x,t) [m]')
% % s=surf(X,T,f);
% % grid
% % xlabel('x [m]')
% % ylabel('t [s]')
% % zlabel('f(x,t) [m]')
% % zlim([min(f(:)),max(f(:))])
% % s.EdgeColor = 'none';
% 
% subplot(2,2,3)
% h = plot(x,f(1,:));
% grid
% xlabel('x [m]')
% ylabel('f(x,t) [m]')
% ht = title('t=0 s');
% ylim([min(f(:)),max(f(:))])
% for i=2:length(t)
%     pause(.01)
%     if ~ishandle(h)
%          break % Arrete l'animation si la fenetre est fermee
%      end
%      set(h,'YData',f(i,:))
%      set(ht,'String',sprintf('t=%0.2f s',t(i)))
%  end
% %}