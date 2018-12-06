% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=0;
tFin=100*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             50        % dt
             1e-5      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);



rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

vy0T=omega*rT
vy0L=omega*rL

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',  'Lune',  'Apollo13' };
variables = [ rT        rL       1e10      % x0
              0         0        0         % y0
              0         0        0         % z0
              0         0        1200      % vx0
              vy0T      vy0L     0         % vy0
              0         0        0         % vz0
              mT        mL       1         % m
              RT        RL       1.95      % R
              0         0        0      ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);


%% Simulations %%
%%%%%%%%%%%%%%%%%

cmd = './performance';
system(cmd);


%% Analyse %%
%%%%%%%%%%%%


% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3

data=load('simulations/deuxCorps.out');
t  = data(:,1);
emec= data(:,3);
xT = data(:,5);
yT = data(:,6);
xL = data(:,11);
yL = data(:,12);


%% Figures %%
%%%%%%%%%%%%%

figure
plot(xT,yT,xL,yL)
% angle=linspace(0,2*pi);
% plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
% hold on
% plot(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), 'k+')
% hold on
% plot(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), 'r+')
% hold off
% xlabel('x [m]')
% ylabel('y [m]')
grid on
axis equal
figure
plot(t,emec)
% figure
% plot3(t, xT, yT, 'k', t, xL, yL, 'r', t, xA, yA)
% xlabel('t [s]')
% ylabel('x [m]')
% zlabel('y [m]')
% plot(dt.*dt,error,'k+')
% xlabel('(\Deltat)^2 [s^2]')
% ylabel('\theta(t_{fin}) [rad]')
% set(gca,'fontsize',15);
% title('$\Omega$=1$\omega_0$  $d$=0.04  $\kappa$=0', 'Fontweight','normal','Interpreter','latex');
% grid on
% print('figures/etudeConvDt', '-depsc');

clear all;