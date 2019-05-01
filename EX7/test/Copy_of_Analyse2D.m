%% Donnees %%
%%%%%%%%%%%%%

%%% Chargement des donnees %%%
fichier = 'output';
data = load([fichier,'_u.out']);
x = data(:,1);
data = load([fichier,'_f.out']);
t=data(:,1);
f_ = data(:,2:end); % amplitudes en lignes 

%%% determintaion des parametres numerique %%%
Nx=length(x);
Ny=size(f_,2)/Nx;
y = min(x(2:Nx)-x(1:Nx-1))*(0:Ny-1)'; % pour reconstruire les positions y

%%% passage de la represenation en ligne a celle en grille 2D %%%
[X,Y]=meshgrid(x,y);
f=zeros(Ny,Nx,length(t));
for i=1:length(t)
    f(:,:,i)=reshape(f_(i,:),Ny,Nx);
end

%% Figures %%
%%%%%%%%%%%%%

i=800;

%%% plot %%%
fig1=figure('Name',['Analyse de ' fichier]);
h = surf(X,Y,f(:,:,i));

%%% legendes %%%
xlabel('x [m]')
ylabel('y [m]')
zlabel('f(x,y,t) [m]')
ht = title(sprintf('t=%0.2f s',t(i)));

%%% axes %%%
Axy=min(x(end)-x(1),y(end)-y(1));
Az=max(f(:))-min(f(:));
set(gca,'DataAspectRatio',[Axy Axy Az*2])
xlim([x(1),x(end)])
ylim([y(1),y(end)])
zlim([min(f(:)),max(f(:))])

%%% options supplementaires %%%
% colorbar                            % pour avoir la barre de couleur (ralentit l'animation)
caxis([min(f(:)),max(f(:))])        % une couleur = une valeur donnee pour toute l'animation
% h.EdgeColor = 'none';               % enleve les lignes 
% h.FaceAlpha=0.8;                    % rend semi-transparent
% h.FaceLighting='gouraud';           % eclairage different
h.FaceColor='interp';               % interpole les couleurs entre les points de maillages

print(fig1,sprintf('vague_2d_%d',i), '-djpeg');

%%% animation %%%
% w=waitforbuttonpress;
% for i=2:length(t)
%     pause(.01)
%     if ~ishandle(h)
%         break % Arrete l'animation si la fenetre est fermee
%     end
%     set(h,'ZData',f(:,:,i))
%     set(ht,'String',sprintf('t=%0.2f s',t(i)))
% end

    