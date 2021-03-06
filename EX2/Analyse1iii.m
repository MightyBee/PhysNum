% Nom du fichier d'output a analyser
repertoireOut = 'stabilite';
schema = 'Euler';
repertoireExe = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice2'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration.in'; % Nom du fichier d'entree de base


nsteps = [700,700,700]; % Nombre d'iterations entier de 10^2 a 10^4
nsimul = size(nsteps,2); % Nombre de simulations a faire
tfin = 1.09e-7; % TODO: Remplacer la valeur de tfin
dt = tfin ./ nsteps;

paramstr = 'nsteps'; % Nom du parametre a scanner
param = nsteps; % Valeurs du parametre a scanner

%% Simulations %%
%%%%%%%%%%%%%%%%%
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
  output{i} = [repertoireOut,'/',schema,'/',paramstr, '=', num2str(param(i)), '.out'];
  % Execution du programme en lui envoyant la valeur a scanner en argument
  cmd = sprintf('%s%s %s schema=%s %s=%.15g output=%s', repertoireExe, executable, input,schema, paramstr, param(i), output{i});
  disp(cmd)
  system(cmd);
end


output1 = load(output{1});
output2 = load(output{2});
output3 = load(output{3});

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

t3 = output3(:,1);
x3 = output3(:,2);
y3 = output3(:,3);
vx3 = output3(:,4);
vy3 = output3(:,5);
energy3 = output3(:,6);

clear output1
clear output2
clear output3

v0=4e5;
omega=(3*1.6022/1.6726)*1e8;
x0=-1.39e-3;%-v0/omega;


%t = t3;
%x_th = x3;
%y_th = y3;
%
%for i = 1:size(t3)
%  x_th(i) = v0/omega*(1-cos(omega*t(i)))+x0; % TODO: Entrer la vraie solution analytique en fonction du temps
%  y_th(i) = v0/omega*sin(omega*t(i)); % TODO: Entrer la vraie solution analytique en fonction du temps
%end
t=linspace(0,tfin,1000);
x_th = v0/omega*(1-cos(omega*t))+x0; % TODO: Entrer la vraie solution analytique en fonction du temps
y_th = v0/omega*sin(omega*t); % TODO: Entrer la vraie solution analytique en fonction du temps
vx_th=v0*sin(omega*t);
vy_th=v0*cos(omega*t);


%r1=x1;
%for i = 1:size(t1,1)
%  r1(i)=sqrt(x1(i)*x1(i)+y1(i)*y1(i));
%end
% Figures

fig1=figure('Position',[20,20,800,400])
subplot(1,2,1);
plot(x_th,y_th, '--k',x1,y1);%,x2,y2,x3,y3)
axis equal
grid on
xlabel('x [m]','Fontsize', 14)
ylabel('y [m]','Fontsize', 14)
set(gca,'fontsize',14);

subplot(1,2,2);
plot(vx_th,vy_th,'--k',vx1,vy1,'r');%,vx2,vy2,vx3,vy3)
axis equal
grid on
xlabel('v_x [m/s]','Fontsize', 14)
ylabel('v_y [m/s]','Fontsize', 14)
%lgd1=legend('Solution analytique',['Euler-Cromer : ',num2str(nsteps(1)),' steps'],['Euler-Cromer : ',num2str(nsteps(2)),' steps'],['Euler-Cromer : ',num2str(nsteps(3)),' steps'],'Location','northwest')
set(gca,'fontsize',14);
%set(lgd1,'FontSize', 10);
print(fig1,['etudeStabiliteXY',schema], '-depsc');

fig3=figure('Position',[20,20,700,400])
plot(t,x_th,'--',t1,x1,'-',t2,x2,'-',t3,x3,'-')
grid on
xlabel('t [s]','Fontsize', 13)
ylabel('x(t) [m]','Fontsize', 13)
lgd3=legend('Solution analytique',[num2str(nsteps(1)),' steps'],[num2str(nsteps(2)),' steps'],[num2str(nsteps(3)),' steps'],'Location','northeast')
set(gca,'fontsize',13);
set(lgd3,'FontSize', 11);

%print(fig3,['etudeStabilitePos',schema], '-depsc');

fig4=figure('Position',[20,20,700,400])
semilogy(t1,energy1,t2,energy2,t3,energy3)
grid on
xlabel('t [s]','Fontsize', 13)
ylabel('E [J]','Fontsize', 13)
lgd4=legend([num2str(nsteps(1)),' steps'],[num2str(nsteps(2)),' steps'],[num2str(nsteps(3)),' steps'],'Location','northwest')
set(lgd4,'FontSize', 12);
set(gca,'fontsize',13);

%print(fig4,['etudeStabiliteEnergy',schema], '-depsc');

%print([repertoireOut,'/',filename], '-depsc');
