% Nom du fichier d'output a analyser
repertoireOut = 'application1';
filename = 'RK2';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile1 = [repertoireOut,'/',filename, '.out'];
cmd = sprintf('%s%s %s E=6e4 output=%s', repertoireExe, executable, input, outputFile1);
disp(cmd)
system(cmd);
outputFile2 = [repertoireOut,'/',filename, 'E0.out'];
cmd = sprintf('%s%s %s E=0 output=%s', repertoireExe, executable, input, outputFile2);
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
fig1=figure('Position',[20,20,800,400])
subplot(1,2,1)
plot(x1,y1,x2,y2,'k--')
axis equal
grid on
xlabel('x [m]', 'fontsize', 13)
ylabel('y [m]', 'fontsize', 13)
set(gca,'fontsize',13)
lgd = legend('E=6E4 V/M', 'E=0', 'Location', 'southeast');
set(lgd,'FontSize', 13);

subplot(1,2,2)
plot(vx1,vy1,vx2,vy2,'k--')
axis equal
grid on
xlabel('v_x [m/s]', 'fontsize', 13)
ylabel('v_y [m/s]', 'fontsize', 13)
set(gca,'fontsize',13)

print(fig1,[repertoireOut,'/',filename,'XY1'], '-depsc');


fig2=figure('Position',[20,20,800,400])
plot(t1,energy1)
axis auto
grid on
xlabel('t [s]', 'fontsize', 15)
ylabel('E [J]', 'fontsize', 15)
set(gca,'fontsize',15)
print(fig2,[repertoireOut,'/',filename,'Energy1'], '-depsc');
