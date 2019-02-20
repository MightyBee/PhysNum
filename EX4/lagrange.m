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
             1      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);


rTL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rTL/(mT+mL);
rT=0-rG;
rL=rTL-rG;

dev=1000000
if dev>0
    strDev='_dev';
else
    strDev='';
end
Z=0;%100000

rLag=zeros(5,2);
vLag=zeros(5,2);


rLag(4,1)=rT+rTL/2;
rLag(4,2)=sqrt(3)*rTL/2;
rL4=sqrt(rLag(4,1)^2+rLag(4,2)^2);

rLag(5,1)=rT+rTL/2;
rLag(5,2)=-sqrt(3)*rTL/2;
rL5=sqrt(rLag(5,1)^2+rLag(5,2)^2);

rLag(1,1)=point_lagrange( 321000000, 323000000,1e-6);
rLag(2,1)=point_lagrange( 444000000, 445000000,8e-7);
rLag(3,1)=point_lagrange(-420000000,-360000000,1e-7);

omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);
vLag(1,2)=rLag(1,1)*omega;
vLag(2,2)=rLag(2,1)*omega;
vLag(3,2)=rLag(3,1)*omega;
v0L4=rL4*omega;
vLag(4,1)=-v0L4*rLag(4,2)/rL4;
vLag(4,2)=v0L4*rLag(4,1)/rL4;
v0L5=rL5*omega;
vLag(5,1)=-v0L5*rLag(5,2)/rL5;
vLag(5,2)=v0L5*rLag(5,1)/rL5;

vy0T=omega*rT;
vy0L=omega*rL;

RT=6378100;
RL=1737000;

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',  'Lune', 'L1',          'L2',           'L3',         'L4',                'L5'                  };
variables = [ rT        rL      rLag(1,1)-dev  rLag(2,1)+dev  rLag(3,1)+dev  rLag(4,1)+sqrt(dev)  rLag(5,1)+sqrt(dev)  % x0
              0         0       0              0              0              rLag(4,2)+sqrt(dev)  rLag(5,2)-sqrt(dev)  % y0
              0         0       Z              Z              Z              Z                    Z                    % z0
              0         0       0              0              0              vLag(4,1)            vLag(5,1)            % vx0
              vy0T      vy0L    vLag(1,2)      vLag(2,2)      vLag(3,2)      vLag(4,2)            vLag(5,2)            % vy0
              0         0       0              0              0              0                    0                    % vz0
              mT        mL      10             10             10             10                   10                   % m
              RT        RL      1.95           1.95           1.95           1.95                 1.95                 % R
              0         0       0              0              0              0                    0      ];      % Cx

% iLag=[1 2 3 4 5];
% T1=table(variables(:,1),variables(:,2),variables(:,3),variables(:,4),variables(:,5),variables(:,6),variables(:,7),'VariableNames',varNames,'RowNames',rowNames);

% iLag=[3 4 5];
% T1=table(variables(:,1),variables(:,2),variables(:,iLag(1)+2),variables(:,iLag(2)+2),variables(:,iLag(3)+2),'VariableNames',varNames([1 2 iLag+2]),'RowNames',rowNames);
% nbCorps=size(iLag,2)+2;

iLag=[4];
T1=table(variables(:,1),variables(:,2),variables(:,iLag(1)+2),'VariableNames',varNames([1 2 iLag+2]),'RowNames',rowNames);
nbCorps=size(iLag,2)+2;

config(T0,T1);
change_config(0,'nbCorps',nbCorps);

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
Pot=-G*(mT./r13+mL./r23)-0.5*omega^2*r3.^2;
for i=1:n
    for k=1:n
        if Pot(i,k) < -1.7e6
            Pot(i,k)=nan;
        end
    end
end
iT= r13 > 1e7;
angle=linspace(0,2*pi,10000);

fig1=figure('Position',[50,50,550,400]);
plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
hold on
contourf(x,y,Pot,13)
hold on
plot(rLag(1,1),0,'rp',rLag(2,1),0,'rp',rLag(3,1),0,'rp',rLag(4,1),rLag(4,2),'rp',rLag(5,1),rLag(5,2),'rp','MarkerSize',10)
hold off
xlabel("x' [m]")
ylabel("y' [m]")
cbar = colorbar;
cbar.Label.String = 'Potentiel [m^2s^{-2}]';
set(gca,'fontsize',15);
axis equal
print(fig1,'figures/lagrange/lagrange_Epot', '-depsc');

colorspec = {[0 0.447 0.741]; [0.85 0.325 0.098]; [0.929 0.694 0.125]; [0.494 0.184 0.556]; [0.466 0.674 0.188]; [0.301 0.745 0.933]; [0.635 0.078 0.184]};

% if Z==0
    fig2=figure('Position',[50,50,550,400]);
    if nbCorps>3
        angle=linspace(0,2*pi,10000);
        plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'Color',colorspec{1})
        hold on
        plot(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), '+','Color',colorspec{1})
        hold on
        plot(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), '+','Color',colorspec{2})
        hold on
    end
    for i=1:nbCorps-2
        plot(rLag(iLag(i),1),rLag(iLag(i),2),'p','Color',colorspec{iLag(i)+2})
        hold on
        plot(xLag(:,i).*cos(omega*t)+yLag(:,i).*sin(omega*t), -xLag(:,i).*sin(omega*t)+yLag(:,i).*cos(omega*t), '-','Color',colorspec{iLag(i)+2})
        hold on
    end
    hold off
    xlabel("x' [m]")
    ylabel("y' [m]")
    set(gca,'fontsize',15);
    grid on
    axis equal
    print(fig2,sprintf('figures/lagrange/lagrange_%dpoints%s_trajectoireRprime',nbCorps-2,strDev), '-depsc');
% else
%     figure
% %     angle=linspace(0,2*pi,10000);
% %     plot(rT+RT*cos(angle),RT*sin(angle),'r',rL+RL*cos(angle),RL*sin(angle),'r')
% %     hold on
%     plot3(xT.*cos(omega*t)+yT.*sin(omega*t), -xT.*sin(omega*t)+yT.*cos(omega*t), zT, '+')
%     hold on
%     plot3(xL.*cos(omega*t)+yL.*sin(omega*t), -xL.*sin(omega*t)+yL.*cos(omega*t), zL, '+')
%     hold on
%     for i=1:nbCorps-2
%         plot3(xLag(1,i),yLag(1,i),zLag(1,i),'bp',xLag(:,i).*cos(omega*t)+yLag(:,i).*sin(omega*t), -xLag(:,i).*sin(omega*t)+yLag(:,i).*cos(omega*t), zLag(:,i), 'b-')
%         hold on
%     end
%     hold off
%     xlabel("x' [m]")
%     ylabel("y' [m]")
%     zlabel("z' [m]")
%     grid on
%     axis equal
% end
    

fig3=figure('Position',[50,50,550,400]);
plot(xT,yT,'Color',colorspec{1})
hold on 
plot(xL,yL,'Color',colorspec{2})
hold on
for i=1:nbCorps-2
    plot(xLag(:,i),yLag(:,i),'Color',colorspec{iLag(i)+2})
    hold on
end
hold off
xlabel('x [m]')
ylabel('y [m]')
set(gca,'fontsize',15);
grid on
axis equal
print(fig3,sprintf('figures/lagrange/lagrange_%dpoints%s_trajectoireR',nbCorps-2,strDev), '-depsc');

fig4=figure('Position',[50,50,600,400]);
if nbCorps>3
    plot(t,sqrt(xT.^2+yT.^2+zT.^2),'Color',colorspec{1})
    hold on
    plot(t,sqrt(xL.^2+yL.^2+zL.^2),'Color',colorspec{2})
    hold on
end
for i=1:nbCorps-2
    plot(t,sqrt(xLag(:,i).^2+yLag(:,i).^2+zLag.^2),'Color',colorspec{iLag(i)+2})
    hold on
end
hold off
xlabel('t [s]')
ylabel('r [m]')
set(gca,'fontsize',15);
grid on
print(fig4,sprintf('figures/lagrange/lagrange_%dpoints%s_distance',nbCorps-2,strDev), '-depsc');

