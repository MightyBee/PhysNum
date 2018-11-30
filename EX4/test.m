cd %% ConfigFile %%
%%%%%%%%%%%%%%%%

v=1200;
alpha=pi-10.86*pi/180
vx0=v*cos(alpha);
vy0=v*sin(alpha);

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [0               384748000       314159000       % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             5.972e24        7.3477e8        5809            % m
             6378100         1737000         1.95            % R
             0               0               0            ]; % Cx

T=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T);

change_config(0,'tFin',100*24*60*60);
change_config(0,'dt',60);

cmd = './performance';
% system(cmd);
output = load('simulations/h.out');
% Extraction des quantites d'interet
t = output(:,1);
x1 = output(:,2);
y1 = output(:,3);
% x2 = output(:,8);
% y2 = output(:,9);
x3 = output(:,14);
y3 = output(:,15);
clear output

figure
plot(x1,y1,'+', x3,y3)
grid on 
axis equal
% figure
% plot3(t,x1,y1,'+',t,t, x3,y3)
% grid on
% figure
% plot3(t,x3-x3,y3-y3,'+',t,-x3,-y3)
% grid on
% figure
% plot(t/(24*60*60),sqrt((-x3).^2+(-y3).^2))
% grid on


