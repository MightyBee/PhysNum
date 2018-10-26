% Nom du fichier d'output a analyser
repertoireOut = '';
filename = 'a';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = [repertoireOut,filename, '.out'];
cmd = sprintf('%s%s %s output=%s', repertoireExe, executable, input, outputFile);
disp(cmd)
system(cmd);

% Chargement des donnees
output = load(outputFile);

% Extraction des quantites d'interet
t = output(:,1);
theta = output(:,2);
thetadot = output(:,3);
energy= output(:,4);
P = output(:,5);
x = 0.1*sin(theta)
y= -0.1*cos(theta)
clear output

% Figures

figure
%subplot(2,3,1)
plot(x,y);
%axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')

figure
%subplot(2,3,1)
plot(t,theta);
%axis equal
grid on
xlabel('t [s]')
ylabel('theta [rad]')

figure
%subplot(2,3,2)
plot(thetadot,t)
%axis equal
grid on
xlabel('t [s]')
ylabel('thetadot [rad/s]')

figure
%subplot(2,3,3)
plot(t,x,t,y)
%grid on
%label('tx [s]')
%ylabel('x,y [m]')
%legend('x','y')

%subplot(2,3,4)
%plot(t,vx,t,vy)
%grid on
%xlabel('t [s]')
%ylabel('v_x,v_y [m/s]')
%legend('v_x','v_y')

%subplot(2,3,5)
%plot(t,energy)
%grid on
%xlabel('t [s]')
%ylabel('E [J]')

%subplot(2,3,6)
%plot(t,mu)
%grid on
%xlabel('t [s]')
%ylabel('\mu [J/T]')

%print([repertoireOut,'/',filename], '-depsc');
