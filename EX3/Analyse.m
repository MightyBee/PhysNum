% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'a';
repertoireExe = './'; % Chemin d'acces au code compile
executable = 'Exercice3'; % Nom de l'executable
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
x =  0.1*sin(theta);
y = -0.1*cos(theta);
yAbs = y + 0.03*sin(9.9045*t);
clear output


% Figures

figure
plot3(theta.*cos(t),theta.*sin(t),thetadot);
grid on 


figure
%subplot(2,3,1)
plot(x,y)
axis equal
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
plot(t,thetadot)
%axis equal
grid on
xlabel('t [s]')
ylabel('thetadot [rad/s]')

figure
%subplot(2,3,3)
Esim=plot(t,energy,'b',[0 20],[0.1*9.81*0.1*(1-cos(0.000001)) 0.1*9.81*0.1*(1-cos(0.000001))], 'r')


grid on
xlabel('t [s]')
ylabel('energy[J]')
%legend('x','y')


