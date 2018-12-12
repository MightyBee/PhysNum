% cd %% ConfigFile %%
% %%%%%%%%%%%%%%%%
% 

G=6.674e-11;
rho0=1.2;
tFin=1000*24*3600;

rTL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rTL/(mT+mL);
rT=0-rG;
rL=rTL-rG;
deviation=1000000;
x0L4=rT+rTL/2+sqrt(deviation);
y0L4=sqrt(3)*rTL/2+sqrt(deviation);
rL4=sqrt(x0L4^2+y0L4^2);
x0L5=rT+rTL/2+sqrt(deviation);
y0L5=-sqrt(3)*rTL/2-sqrt(deviation);
rL5=sqrt(x0L5^2+y0L5^2);
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);
v0L4=rL4*omega;
vx0L4=-v0L4*y0L4/rL4;
vy0L4=v0L4*x0L4/rL4;
v0L5=rL5*omega;
vx0L5=-v0L5*y0L5/rL5;
vy0L5=v0L5*x0L5/rL5;

vy0T=omega*rT
vy0L=omega*rL


RT=6378100;
RL=1737000;


% v=1200;
% alpha=pi-10.86*pi/180
% vx0=v*cos(alpha);
% vy0=v*sin(alpha);
% 
% rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
% varNames  = {'Terre',        'Lune',         'Apollo13'        };
% variables = [0               384748000       314159000       % x0
%              0               0               0               % y0
%              0               0               0               % z0
%              0               0               vx0             % vx0
%              0               0               vy0             % vy0
%              0               0               0               % vz0
%              5.972e24        7.3477e8        5809            % m
%              6378100         1737000         1.95            % R
%              0              0               0.3            ]; % Cx
% 
% T=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);
% 
% config(T);
% 
% change_config(0,'tFin',100*24*60*60);
% change_config(0,'dt',60);
% 
% cmd = './performance';
% % system(cmd);
% output = load('simulations/h.out');
% % Extraction des quantites d'interet
% t = output(:,1);
% x1 = output(:,2);
% y1 = output(:,3);
% % x2 = output(:,8);
% % y2 = output(:,9);
% x3 = output(:,14);
% y3 = output(:,15);
% clear output
% 
% figure
% plot(x1,y1,'+', x3,y3)
% grid on 
% axis equal
% % figure
% % plot3(t,x1,y1,'+',t,t, x3,y3)
% % grid on
% % figure
% % plot3(t,x3-x3,y3-y3,'+',t,-x3,-y3)
% % grid on
% % figure
% % plot(t/(24*60*60),sqrt((-x3).^2+(-y3).^2))
% % grid on

L1=point_lagrange( 300000000, 350000000,1);
L2=point_lagrange( 400000000, 500000000,1);
L3=point_lagrange(-420000000,-360000000,1);

n=1000;
[theta,rayon] = meshgrid(linspace(0,2*pi,n),linspace(2e8,5e8,n));
[x,y] = pol2cart(theta,rayon);
r13=sqrt((x-rT).^2+(y).^2);
r23=sqrt((x-rL).^2+(y).^2);
r3 =sqrt(x.^2+y.^2);
Epot=-G*(mT./r13+mL./r23)-0.5*omega^2*r3.^2;
for i=1:n
    for k=1:n
%         if  r23(i,k) < 3.5e7
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
plot(L1,0,'rp',L2,0,'rp',L3,0,'rp','MarkerSize',10)
hold off

x=[linspace(-420000000,-360000000,100000) linspace(321000000, 323000000,100000) linspace(444000000, 445000000,100000)];
d=384748000;
mT=5.972e24;
mL=7.3477e22;
alpha=mL/(mT+mL);
beta=mT/(mT+mL);
rG=alpha*d;
xT=0-rG;
xL=d-rG;
y=d^3*(beta*(x-xT)./abs(x-xT).^3+alpha*(x-xL)./abs(x-xL).^3)-x;
figure
plot(x,y)
grid on
