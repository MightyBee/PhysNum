repertoireOut = 'simulations/';
filename = 'a';
repertoireExe = './'; % Chemin d'acces au code compile
executable = 'Exercice3'; % Nom de l'executable
input = 'configuration.in'; % Nom du fichier d'entree de base
nsimul = 3; % Nombre de simulations a faire


It=[0 20];
dt = [0.008,0.08,0.2];

paramstr = 'dt'; % Nom du parametre a scanner (changer ici 'dt' ou 'Omega' ou autre)
param = dt; % Valeurs du parametre a scanner (changer ici dt ou Omega ou autre)

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = ['simulations/',paramstr, '=', num2str(param(i)), '.out'];
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %s=%.15g output=%s', repertoire, executable, input, paramstr, param(i), output{i});
    disp(cmd)
    system(cmd);
end


output1 = load(output{1});
output2 = load(output{2});
output3 = load(output{3});


theta1 = output1(:,2);
thetadot1 = output1(:,3);
t1 = output1(:,1);
energy1= output1(:,4);
P1 = output1(:,5);

theta2 = output2(:,2);
thetadot2 = output2(:,3);
t2 = output2(:,1);
energy2= output2(:,4);
P2 = output2(:,5);


theta3 = output3(:,2);
thetadot3 = output3(:,3);
t3 = output3(:,1);
energy3= output3(:,4);
P3 = output3(:,5);


clear output1
clear output2
clear output3


%figure
%plot3(theta.*cos(t),theta.*sin(t),thetadot);
%grid on 


figure
%subplot(2,3,1)
plot(theta1,thetadot1,'b',theta2,thetadot2,'r',theta3,thetadot3,'g',0.000001*sin((sqrt(9.81/0.1)*t1)+3.14/2),0.000001*sqrt(9.81/0.1)*cos((sqrt(9.81/0.1)*t1)+3.14/2),'k')
title('theta en fonction de thetadot')
legend('dt=0.008','dt=0.08','dt=0.2','solution analytique')
%axis equal
grid on
xlabel('theta [rad]')
ylabel('thetadot [rad/s]')

figure
%subplot(2,3,1)
plot(t1,theta1,'b',t2,theta2,'r',t3,theta3,'g',t1,0.000001*sin((sqrt(9.81/0.1)*t1)+3.14/2),'k')
%axis equal
title('theta en fonction de t')
legend('dt=0.008','dt=0.08','dt=0.2','solution analytique')
grid on
xlabel('t [s]')
ylabel('theta [rad]')

figure
%subplot(2,3,2)
plot(t1,thetadot1,'b',t2,thetadot2,'r',t3,thetadot3,'g',t1,0.000001*sqrt(9.81/0.1)*cos((sqrt(9.81/0.1)*t1)+3.14/2),'k')
%axis equal
title('thetadot en fonction de t')
legend('dt=0.008','dt=0.08','dt=0.2','solution analytique')
grid on
xlabel('t [s]')
ylabel('thetadot [rad/s]')

figure
plot(t1,energy1,'b',t2,energy2,'r',t3,energy3,'g',[0 20],[0.1*9.81*0.1*(1-cos(0.000001)) 0.1*9.81*0.1*(1-cos(0.000001))], 'k')
grid on
title('energie en fonction de t')
legend('dt=0.008','dt=0.08','dt=0.2','solution analytique','Location','southeast')
xlabel('t [s]')
ylabel('energy[J]')
