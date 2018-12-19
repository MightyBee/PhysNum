% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=0;
tFin=2*24*3600;
output="a.out"

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             1         % dt
             1e-6      % precision
             "false"    % adaptatif
             output    % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);



v0=1200;
r0=314159000;
mT=5.972e24;
h=10000;
RT=6378100;
vMax_th=sqrt(v0^2+2*G*mT*(1/(h+RT)-1/r0));
alpha = pi-asin(vMax_th*(h+RT)/(v0*r0))
vx0=v0*cos(alpha);
vy0=v0*sin(alpha);

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [0               0               r0              % x0
             0               384748000       0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             mT              7.3477e0        5809            % m
             RT              1737000         1.95            % R
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);

%% simulation %%

cmd = "./performance";
disp(cmd)
system(cmd);


%% Analyse %%
%%%%%%%%%%%%%


% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2  


data = load("simulations/"+output); % Chargement du fichier de sortie de la i-ieme simulation
t = data(:,1);
xA = data(:,17);
yA = data(:,18);
xT = data(:,5);
yT = data(:,6);
vxA= data(:,20);
vyA= data(:,21);

v = sqrt(vxA.^2 + vyA.^2);
[value, indice]=min(sqrt(xA.^2+yA.^2));

%% Figures %%
%%%%%%%%%%%%%

colorspec = {[0 0.447 0.741]; [0.85 0.325 0.098]};

fig1=figure('Position',[50,50,400,420]);
% g = patch('Vertices', [xA(:), yA(:); nan nan ], 'Faces', (1:length(xA)+1).', 'FaceVertexCData', [v(:); nan], 'EdgeColor', 'interp', 'Marker','.','MarkerSize',4);
% cbar = colorbar;
% cbar.Label.String = 'v [m/s]';
plot(xA,yA,'Color',colorspec{1})
hold on
plot(xA(indice),yA(indice),'x','Color',colorspec{1})
hold on
plot(RT*cos(linspace(0,2*pi,100000)),RT*sin(linspace(0,2*pi,100000)),'Color',colorspec{2})
hold off
% ylim([-16e7 4e7])
% xlim([-0.2e8 3.3e8])
xlabel('x [m]')
ylabel('y [m]')
set(gca,'fontsize',14);
axis equal
grid on
print(fig1,'figures/vitesssssseee','-depsc');