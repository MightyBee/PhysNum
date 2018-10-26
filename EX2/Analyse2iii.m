% Nom du fichier d'output a analyser
repertoireOut = 'application1';
filename = 'RK2';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = [repertoireOut,'/',filename, '.out'];
cmd = sprintf('%s%s %s E=6e4 output=%s', repertoireExe, executable, input, outputFile);
disp(cmd)
system(cmd);

% Chargement des donnees
output = load(outputFile);

% Extraction des quantites d'interet
t = output(:,1);
x = output(:,2);
y = output(:,3);
vx = output(:,4);
vy = output(:,5);
energy = output(:,6);
mu = output(:,7);

clear output
E=6e4;
B=3;
vE=E/B;
for i = 1:size(t)
  x(i)=x(i)-vE*t(i);
end
% Figures

figure('Position',[20,20,600,400])
plot(x,y)
axis equal
grid on
xlabel('x [m]','fontsize', 13)
ylabel('y [m]','fontsize', 13)
set(gca,'fontsize',13)


print([repertoireOut,'/',filename,'XYvE1'], '-depsc');
