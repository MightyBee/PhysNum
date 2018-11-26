% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'poincare';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base

is=100;
ie=10000;
sampling=200;
%           1    2   3 4 5 6   7    8 9   10   11  12    13  14
theta0=    [0    0   0 0 0 0.5 1e-3 2 0   1e-6 3   2     2   1  3  pi+0.01 pi+0.02 pi+1 pi-1 1];
thetadot0= [1e-2 0.5 2 3 4 4   50   0 150 0    -15 -17.8 -34 0  0  0       0       3    0    1.1];
% taille = 30;
% % theta0=logspace(-2,0.477,taille);
% theta0=[linspace(0,1.35,round(taille/4)) linspace(1.4,1.59,taille-round(taille/4))];
% thetadot0=ones(taille)*1e-3;

w0=sqrt(9.81/0.1); %%% 15.1197 pour theta0=0, thetadot0=1e-1 et d=0.04
Omega=6*w0;
d=0.03;
kappa=0;
dt=2*pi/(sampling*Omega);
tFin=dt*ie*sampling;

aAnalyser = [16 17]
% aAnalyser = 1:taille;
nsimul = size(aAnalyser,2);

%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = aAnalyser
    outputFile{i} = [repertoireOut,filename,'_theta0=',num2str(theta0(i)),'_thetadot0=',num2str(thetadot0(i)),'_i=', num2str(i),'.out'];
    cmd = sprintf('%s%s %s tFin=%.15g d=%.15g Omega=%.20g kappa=%.15g theta0=%.15g thetadot0=%.15g dt=%.15g output=%s sampling=%d', repertoireExe, executable, input, tFin, d, Omega, kappa, theta0(i), thetadot0(i), dt, outputFile{i}, sampling);
    disp(cmd)
    system(cmd);
end


figure('Position',[50,50,600,400]);


for i= aAnalyser
    % Chargement des donnees
    output = load(outputFile{i});

    % Extraction des quantites d'interet
    theta = output(:,2);
    thetadot = output(:,3);
    clear output
    
    % Figures
    plot(mod(theta(is:ie)+pi,2*pi)-pi, thetadot(is:ie),'.','Markersize',4);
    hold on;
end
hold off;
xlabel('\theta [rad]')
ylabel('\omega [rad/s]')
% xlim([-pi pi])
% ylim([-35,45])
set(gca,'fontsize',15);
lgd = legend('$\theta_0 = \pi + 0.01$ rad', '$\theta_0=\pi + 0.02$ rad' , 'Location', 'northeast','Interpreter','latex');
set(lgd,'FontSize', 13);
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
% lgd=legend('10^{-2}','$\Dot{\theta} =  10^{-2}$','$\Dot{\theta} =  10^{-2}$','$\Dot{\theta} =  10^{-2}$','$\Dot{\theta} =  10^{-2}$');
% set(lgd,'fontsize',14,'Location','northwest');
grid on
print('figures/poincare', '-depsc');



















% % Chargement des donnees
% output = load(outputFile{1});
% 
% % Extraction des quantites d'interet
% t=output(100:10000,1);
% theta = output(100:10000,2);
% thetadot = output(100:10000,3);
% clear output
% figure
% plot3(theta.*cos(t/Omega),theta.*sin(t/Omega),thetadot);
% grid on 

% figure
% for i= aAnalyser
%     % Chargement des donnees
%     output = load(outputFile{i});
% 
%     % Extraction des quantites d'interet
%     theta = output(:,2);
%     thetadot = output(:,3);
%     clear output
%     
%     % Figures
%     plot(0.1*sin(theta(is:ie)), 0.1*thetadot(is:ie).*sin(theta(is:ie)),[colors(i),'.']);
%     hold on;
%     plot(-0.1*cos(theta(is:ie)), -0.1*thetadot(is:ie).*cos(theta(is:ie)),[colors(i),'.']);
%     hold on;
% end
% hold off;
% xlabel('position [m]')
% ylabel('vitesse [m/s]')
% lgd = legend('coordonnée x', 'coordonné y', 'Location', 'southeast');
% set(lgd,'FontSize', 13);
% set(gca,'fontsize',15);
% grid on
