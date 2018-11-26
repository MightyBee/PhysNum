% Nom du fichier d'output a analyser
repertoireOut = 'simulations/';
filename = 'poincare';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice3'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base


sampling=100;
theta0=3;
thetadot0=0;
w0=sqrt(9.81/0.1)
Omega=2*w0; %%% 15.1197 pour theta0=0, thetadot0=1e-1 et d=0.04
d=0.05;
kappa=0.1;
dt=2e-4;
tFin=200*2*pi/Omega;

nsimul = 2;
%% Simulations %%
%%%%%%%%%%%%%%%%%
% Execution du programme en lui envoyant la valeur a scanner en argument
outputFile = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    outputFile{i} = [repertoireOut,filename,'_theta0=',num2str(theta0),'_thetadot0=',num2str(thetadot0),'_', num2str(i), '.out'];
    cmd = sprintf('%s%s %s tFin=%.15g d=%.15g kappa=%.15g Omega=%.15g theta0=%.15g thetadot0=%.15g dt=%.15g output=%s sampling=%d', repertoireExe, executable, input, tFin, d, kappa, Omega, theta0, thetadot0, dt, outputFile{i}, sampling);
    disp(cmd)
    system(cmd);
    theta0=theta0+1e-8;
end



% Chargement des donnees
output = load(outputFile{1});

% Extraction des quantites d'interet
t1 = output(:,1);
theta1 = output(:,2);
thetadot1 = output(:,3);
clear output

output = load(outputFile{2});
% Extraction des quantites d'interet
t2 = output(:,1);
theta2 = output(:,2);
thetadot2 = output(:,3);
clear output


dist=zeros(size(t1,1),1);
for i = 1:size(t1)
    dist(i)=sqrt((thetadot1(i)-thetadot2(i))^2+Omega^2*(theta1(i)-theta2(i))^2);
end

k=round(12/(sampling*dt));
tFit=t1(1:k);
dFit=log(dist(1:k));
X=[ones(k,1) tFit];
coef=X\dFit;
y=X*coef;

% Figures
fig1=figure('Position',[50,50,700,450]);
plot(t1,theta1,'-',t2,theta2,'-');
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
xlabel('t [s]')
ylabel('$\theta$ [rad]','Interpreter','latex')
set(gca,'fontsize',15);
lgd = legend('$\theta_0 = 1$ rad', '$\theta_0=3$ rad' , 'Location', 'northeast','Interpreter','latex');
set(lgd,'FontSize', 13);
grid on
% rectangle('Position',[8 3 9 37]);
% axes('position',[.2 .55 .35 .35])
% box on % put box around new pair of axes
% indexOfInterest = (t1 < 16) & (t1 > 10); % range of t near perturbation
% plot(t1(indexOfInterest),theta1(indexOfInterest),t2(indexOfInterest),theta2(indexOfInterest)) % plot on new axes
% axis tight
% grid on
print(fig1,'figures/lyapunovThetaInstableAmorti', '-depsc');

fig2=figure('Position',[50,50,700,450]);
plot(t1,log(dist),'-');
hold on
semilogy(tFit,y,'-','Linewidth',1.5);
xlabel('t [s]')
ylabel('ln(d) [rad/s]')
set(gca,'fontsize',15);
title(['$\Omega$=' num2str(Omega/w0) '$\omega_0$  d=' num2str(d) '  $\kappa$=' num2str(kappa)], 'Fontweight','normal','Interpreter','latex');
lgd = legend('ln(d(t))',['pente du fit linéaire : ',num2str(coef(2))] , 'Location', 'southeast');
set(lgd,'FontSize', 13);
grid on
print(fig2,'figures/lyapunovOrbiteInstableAmorti', '-depsc');