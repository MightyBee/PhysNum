% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%


rL=384748000;
G=6.674e-11;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

alpha = 0.3 %0.2138; %pi-asin(vMax_th*(h+RT)/(v0A*r0));
v0A=1200;
vx0A=v0A*cos(alpha);
vy0A=v0A*sin(alpha)+omega*rA;

vy0T=omega*rT
vy0L=omega*rL

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',  'Lune',  'Apollo13' };
variables = [ rT        rL       rA        % x0
              0         0        0         % y0
              0         0        0         % z0
              0         0        vx0A      % vx0
              vy0T      vy0L     vy0A      % vy0
              0         0        0         % vz0
              mT        mL       5809      % m
              RT        RL       1.95      % R
              0         0        0      ]; % Cx

T=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T);

change_config(0,'tFin',100*24*60*60);
change_config(0,'precision',0.00001);
change_config(0,'adaptatif','true');


%% Simulations %%
%%%%%%%%%%%%%%%%%

cmd = './performance configuration0.in 1 output=simulations/deuxCorps.out';
system(cmd);


%% Analyse %%
%%%%%%%%%%%%


% 1 2    3  4  5  6   7   8     9  10 11  12  13  14    15 16 17  18  19  20
% t acc  x1 y1 z1 vx1 vy1 vz1   x2 y2 z2  vx2 vy2 vz2   x3 y3 z3  vx3 vy3 vz3

data=load('simulations/deuxCorps.out');
t  = data(:,1);
xT = data(:,3);
yT = data(:,4);
xL = data(:,9);
yL = data(:,10);
xA = data(:,15);
yA = data(:,16);

disp(size(t));


if(min(sqrt((xA-xT).^2+(yA-yT).^2))<=RT || min(sqrt((xA-xL).^2+(yA-yL).^2))<=RL)
    iFin=1;
    while sqrt((xA(iFin)-xT(iFin)).^2+(yA(iFin)-yT(iFin)).^2)>RT && sqrt((xA(iFin)-xL(iFin)).^2+(yA(iFin)-yL(iFin)).^2)>RL
        iFin=iFin+1;
    end
    t  = t(1:iFin);
    xT  = xT(1:iFin);
    yT  = yT(1:iFin);
    xL  = xL(1:iFin);
    yL  = yL(1:iFin);
    xA  = xA(1:iFin);
    yA  = yA(1:iFin);
end

%% Figures %%
%%%%%%%%%%%%%

figure
% plot(xT,yT,xL,yL,xA,yA)
angle=linspace(0,2*pi);
plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
hold on
plot(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), 'k+')
hold on
plot(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), 'r+')
hold on
plot(xA.*cos(omega*t)+yA.*sin(omega*t), -xA.*sin(omega*t)+yA.*cos(omega*t), 'b')
xlabel('x [m]')
ylabel('y [m]')
grid on
axis equal
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