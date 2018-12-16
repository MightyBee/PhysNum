% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%

nbCorps=7;
G=6.674e-11;
tFin=365*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [nbCorps   % nbCorps
             tFin      % tFin
             G         % G
             0         % rho0
             7238.2    % lambda
             50        % dt
             1e-5      % precision
             "true"    % adaptatif
             "lagrange.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);


rTL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rTL/(mT+mL);
rT=0-rG;
rL=rTL-rG;

deviation=10;%1000000
Z=0;%100000

x0L4=rT+rTL/2+sqrt(deviation);
y0L4=sqrt(3)*rTL/2+sqrt(deviation);
rL4=sqrt(x0L4^2+y0L4^2);

x0L5=rT+rTL/2+sqrt(deviation);
y0L5=-sqrt(3)*rTL/2-sqrt(deviation);
rL5=sqrt(x0L5^2+y0L5^2);

L1=point_lagrange( 321000000, 323000000,1e-6)+deviation;
L2=point_lagrange( 444000000, 445000000,8e-7)+deviation;
L3=point_lagrange(-420000000,-360000000,1e-7)-deviation;

omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);
v0L1=L1*omega;
v0L2=L2*omega;
v0L3=L3*omega;
v0L4=rL4*omega;
vx0L4=-v0L4*y0L4/rL4;
vy0L4=v0L4*x0L4/rL4;
v0L5=rL5*omega;
vx0L5=-v0L5*y0L5/rL5;
vy0L5=v0L5*x0L5/rL5;

vy0T=omega*rT;
vy0L=omega*rL;

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',  'Lune',  'L4',  'L5',  'L1',  'L2',  'L3'        };
variables = [ rT        rL       x0L4   x0L5   L1     L2     L3        % x0
              0         0        y0L4   y0L5   0      0      0         % y0
              0         0        Z      Z      Z      Z      Z         % z0
              0         0        vx0L4  vx0L5  0      0      0         % vx0
              vy0T      vy0L     vy0L4  vy0L5  v0L1   v0L2   v0L3      % vy0
              0         0        0      0      0      0      0         % vz0
              mT        mL       10     10     10     10     10        % m
              RT        RL       1.95   1.95   1.95   1.95   1.95      % R
              0         0        0.3    0.3    0.3    0.3    0.3    ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),variables(:,4),variables(:,5),variables(:,6),variables(:,7),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);


%% Simulations %%
%%%%%%%%%%%%%%%%%
if nbCorps==3 && Z==0
    system('./performance');
else
    system('./Exercice4');
end

%% Analyse %%
%%%%%%%%%%%%%

data = load('simulations/lagrange.out');

% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2    x3 y3 z3  vx3 vy3 vz3

t = data(:,1);
Pt = data(:,4);
xT = data(:,5);
yT = data(:,6);
zT = data(:,7);
xL = data(:,11);
yL = data(:,12);
zL = data(:,13);
xLag=zeros(size(t,1),nbCorps-2);
yLag=zeros(size(t,1),nbCorps-2);
zLag=zeros(size(t,1),nbCorps-2);
for i=1:nbCorps-2
    xLag(:,i) = data(:,11+6*i);
    yLag(:,i) = data(:,12+6*i);
    zLag(:,i) = data(:,13+6*i);
end

nsteps=size(t,1)-1;

%% Figures %%
%%%%%%%%%%%%%
% 

n=1000;
[theta,rayon] = meshgrid(linspace(0,2*pi,n),linspace(2e8,5e8,n));
[x,y] = pol2cart(theta,rayon);
r13=sqrt((x-rT).^2+(y).^2);
r23=sqrt((x-rL).^2+(y).^2);
r3 =sqrt(x.^2+y.^2);
Epot=-G*(mT./r13+mL./r23)-0.5*omega^2*r3.^2;
for i=1:n
    for k=1:n
        if Epot(i,k) < -1.7e6
            Epot(i,k)=nan;
        end
    end
end
iT= r13 > 1e7;
angle=linspace(0,2*pi,10000);
figure
plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
hold on
contourf(x,y,Epot,13)
hold on
plot(L1,0,'rp',L2,0,'rp',L3,0,'rp',x0L4,y0L4,'rp',x0L5,y0L5,'rp','MarkerSize',10)
hold off

if Z==0
    figure
    angle=linspace(0,2*pi,10000);
    plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
    hold on
    plot(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), 'k+')
    hold on
    plot(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), 'r+')
    hold on
    for i=1:nbCorps-2
        plot(xLag(1,i),yLag(1,i),'bp',xLag(:,i).*cos(omega*t)+yLag(:,i).*sin(omega*t), -xLag(:,i).*sin(omega*t)+yLag(:,i).*cos(omega*t), 'b-')
        hold on
    end
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    grid on
    axis equal
else
    figure
%     angle=linspace(0,2*pi,10000);
%     plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
%     hold on
    plot3(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), zT, 'k+')
    hold on
    plot3(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), zL, 'r+')
    hold on
    for i=1:nbCorps-2
        plot3(xLag(1,i),yLag(1,i),zLag(1,i),'bp',xLag(:,i).*cos(omega*t)+yLag(:,i).*sin(omega*t), -xLag(:,i).*sin(omega*t)+yLag(:,i).*cos(omega*t), zLag(:,i), 'b-')
        hold on
    end
    hold off
    xlabel('x [m]')
    ylabel('y [m]')
    zlabel('z [m]')
    grid on
    axis equal
end
    

figure
plot(xT,yT,xL,yL)
hold on
for i=1:nbCorps-2
    plot(xLag(:,i),yLag(:,i))
    hold on
end
hold off
xlabel('x [m]')
ylabel('y [m]')
grid on
axis equal

figure
plot(t,sqrt(xT.^2+yT.^2+zT.^2),t,sqrt(xL.^2+yL.^2+zL.^2))
hold on
for i=1:nbCorps-2
    plot(t,sqrt(xLag(:,i).^2+yLag(:,i).^2+zLag.^2))
    hold on
end
hold off
xlabel('t [s]')
ylabel('r [m]')
grid on

