% Nom du fichier d'output a analyser
repertoireOut = 'application2';
filename1 = 'proton';
filename2 = 'antiproton'
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile1 = [repertoireOut,'/simulations/',filename1, '.out'];
outputFile2 = [repertoireOut,'/simulations/',filename2, '.out'];
cmd = sprintf('%s%s %s Kappa=100  q=1.6022e-19 x0=-1.39e-3 output=%s', repertoireExe, executable, input, outputFile1);
disp(cmd)
system(cmd);
cmd = sprintf('%s%s %s Kappa=100  q=-1.6022e-19 x0=1.39e-3 output=%s', repertoireExe, executable, input, outputFile2);
disp(cmd)
system(cmd);

% Chargement des donnees
output1 = load(outputFile1);
output2 = load(outputFile2);

% Extraction des quantites d'interet
t1 = output1(:,1);
x1 = output1(:,2);
y1 = output1(:,3);
vx1 = output1(:,4);
vy1 = output1(:,5);
energy1 = output1(:,6);
mu1 = output1(:,7);

t2 = output2(:,1);
x2 = output2(:,2);
y2 = output2(:,3);
vx2 = output2(:,4);
vy2 = output2(:,5);
energy2 = output2(:,6);
mu2 = output2(:,7);

clear output1
clear output2

% Figures

fig1=figure('Position',[20,20,500,400])
subplot(1,2,1)
plot(x1,y1,'b',x2,y2,'r')
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')

subplot(1,2,2)
plot(vx1,vy1,'b',vx2,vy2,'r')
axis equal
grid on
xlabel('v_x [m/s]')
ylabel('v_y [m/s]')

print(fig1,[repertoireOut,'/comparaisonTrajectoire2'], '-depsc');


fig2=figure('Position',[20,20,400,400])
plot(t1,y1,'b',t2,y2,'r')
grid on
xlabel('t [s]')
ylabel('y [m]')
legend('proton', 'antiproton')

print(fig2,[repertoireOut,'/comparaisonY2'], '-depsc');

fig3=figure('Position',[20,20,800,400])
plot(t1,mu1,'b',t2,mu2,'r')
grid on
xlabel('t [s]')
ylabel('\mu [J/T]')

print(fig3,[repertoireOut,'/comparaisonMu2'], '-depsc');
