%% Chargement des resultats %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = 'output';
data = load([output '_T.out']);
N = sqrt(length(data)); % nb de noeuds
Y = data(1:N,2);
X = data(1:N:N*N,1);
T = reshape(data(:,3),N,N)'; % 1D -> 2Dclos
h = (max(X)-min(X))/(N-1);
data = load([output '_P.out']);
t = data(:,1);
Pc = data(:,2);
Pf = data(:,3);
Ptot = data(:,4);
kappa=1.2;

%% Analyse %%
%%%%%%%%%%%%%
% TODO : calcul du flux de chaleur au centre des cellules du maillage
jx=-kappa.*(T(2:N,:)-T(1:N-1,:))./h; % N-1 x N
jy=-kappa.*(T(:,2:N)-T(:,1:N-1))./h; % N x N-1
Xmid = X(1:N-1)+h/2;
Ymid = Y(1:N-1)+h/2;
jxc = (jx(:,1:N-1)+jx(:,2:N))./2;
jyc = (jy(1:N-1,:)+jy(2:N,:))./2;
jnorm = sqrt(jxc.^2+jyc.^2);



%% Figures %%
%%%%%%%%%%%%%
% Temperature :
figure
contourf(X,Y,T',15,'LineStyle','None'), hold on
stride = 2; % (affiche 1 point sur [stride] dans chaque dimension)
quiver(Xmid(1:stride:end,1:stride:end),Ymid(1:stride:end,1:stride:end),jxc(1:stride:end,1:stride:end)',jyc(1:stride:end,1:stride:end)','k')
xlabel('x [m]')
ylabel('y [m]')
title('T(x,y) [°C]')
colorbar
axis equal

% Flux de chaleur :
figure
contourf(Xmid,Ymid,jnorm',30,'LineStyle','None')
xlabel('x [m]')
ylabel('y [m]')
title('|j|(x,y) [W/m]')
colorbar
axis equal

% Puissance :
figure
plot(t, Pc, t, Pf, t, Pc + Pf, t, Ptot)
xlabel('t [s]')
ylabel('P [W]')
legend('P_c', 'P_f', 'P_c+P_f', 'P_{tot}')
grid on

