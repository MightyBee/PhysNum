rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Soleil',  'Terre',        'Lune',         'Halley'        };
variables = [0          147098074000    147482474000    -5285292747000  % x0
             0          0               0               0               % y0
             0          0               0               0               % z0
             0          0               0               0               % vx0
             0          30287           31282           810             % vy0
             nan          0               0               0               % vz0
             1.9891e30  5.9736e24       7.3477e22       1014e11         % m
             696342000  6371000         1736000         15000           % R
             0          0               0               0            ]; % Cx

T=table(variables(:,1),variables(:,2),variables(:,3),variables(:,4),'VariableNames',varNames,'RowNames',rowNames);

config(T);
% % C = change_config(0,"dt",2000);
% cmd = './performance';
% % system(cmd);
% output = load('simulations/f.out');
% % Extraction des quantites d'interet
% t = output(:,1);
% x1 = output(:,2);
% y1 = output(:,3);
% x2 = output(:,8);
% y2 = output(:,9);
% x3 = output(:,14);
% y3 = output(:,15);
% x4 = output(:,20);
% y4 = output(:,21);
% clear output
% index = (t<76.01*365.25636567*24*60*60);
% 
% figure
% plot(x1(index),y1(index),'+',x2(index),y2(index), x3(index),y3(index), x4(index),y4(index))
% grid on 
% axis equal
% figure
% plot3(t(index),x1(index),y1(index),'+',t(index),x4(index),y4(index),t(index), x3(index),y3(index))
% grid on
% figure
% plot3(t(index),x3(index)-x3(index),y3(index)-y3(index),'+',t(index),x4(index)-x3(index),y4(index)-y3(index))
% grid on
% figure
% plot(t(index)/(24*60*60),sqrt((x4(index)-x3(index)).^2+(y4(index)-y3(index)).^2))
% grid on


% Parametres physiques 
nbCorps   = 3 
tFin      = 3153600000 
G         = 6.674e-11 
rho0      = 0 
lambda    = 100 

% Parametres numeriques 
dt        = 0 
precision = 0.00001 
adaptatif = false 
output    = simulations/h.out 
sampling  = 10 
