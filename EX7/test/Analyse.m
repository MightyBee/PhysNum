%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fichier = 'output';
data = load([fichier,'_u.out']);
x = data(:,1);
u = data(:,2);
data = load([fichier,'_f.out']);
t=data(:,1);
f_ = data(:,2:end);

Nx=length(x);
Ny=size(f_,2)/Nx;
y = min(x(2:Nx)-x(1:Nx-1))*(0:Ny-1)';
[X,Y]=meshgrid(x,y);
f=zeros(Ny,Nx,length(t));
for i=1:length(t)
   f(:,:,i)=reshape(f_(i,:),Ny,Nx); 
end

%% Figures %%
%%%%%%%%%%%%%
figure('Name',['Analyse de ' fichier])
h = surf(X,Y,f(:,:,1));
grid
xlabel('x [m]')
ylabel('y [m]')
zlabel('f(x,t) [m]')
ht = title('t=0 s');
hl = get(gca,'DataAspectRatio');
% if h(3)==1
%     set(gca,'DataAspectRatio',[1 1 1/max(hl(1:2))])
% else
    set(gca,'DataAspectRatio',[1 1 hl(3)])
% end
zlim([min(f(:)),max(f(:))])
w=waitforbuttonpress;


for i=2:length(t)
    pause(.01)
    if ~ishandle(h)
        break % Arrete l'animation si la fenetre est fermee
    end
    set(h,'ZData',f(:,:,i))
    set(ht,'String',sprintf('t=%0.2f s',t(i)))
end
% Clear flag and make it invisible.
set(handles.chkFinishNow, 'Value', 0);

    