% Noms des fichiers d'output a analyser
filename1 = '200_steps.out';
filename2 = '400_steps.out';
filename3 = '800_steps.out';
filename4 = '1600_steps.out';
filename5 = '3200_steps.out';
filename6 = '6400_steps.out';


% Chargement des donnees
data1 = load(filename1);
data2 = load(filename2);
data3 = load(filename3);
data4 = load(filename4);
data5 = load(filename5);
data6 = load(filename6);

% Extraction des quantites d'interet
% (Le code c++ ecrit t, z(t) et v(t) en colonnes.)
t1 = data1(:,1);
z1 = data1(:,2);
v1 = data1(:,3);
t2 = data2(:,1);
z2 = data2(:,2);
v2 = data2(:,3);
t3 = data3(:,1);
z3 = data3(:,2);
v3 = data3(:,3);
t4 = data4(:,1);
z4 = data4(:,2);
v4 = data4(:,3);
t5 = data5(:,1);
z5 = data5(:,2);
v5 = data5(:,3);
t6 = data6(:,1);
z6 = data6(:,2);
v6 = data6(:,3);

% Figures
figure('NumberTitle', 'Off', 'Name', [filename1 ': z(t)'])
plot(t1, z1, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename1 ': v(t)'])
plot(t1, v1, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename1 ': v(z)'])
plot(z1, v1, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


figure('NumberTitle', 'Off', 'Name', [filename2 ': z(t)'])
plot(t2, z2, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename2 ': v(t)'])
plot(t2, v2, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename2 ': v(z)'])
plot(z2, v2, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


figure('NumberTitle', 'Off', 'Name', [filename3 ': z(t)'])
plot(t3, z3, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename3 ': v(t)'])
plot(t3, v3, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename3 ': v(z)'])
plot(z3, v3, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


figure('NumberTitle', 'Off', 'Name', [filename4 ': z(t)'])
plot(t4, z4, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename4 ': v(t)'])
plot(t4, v4, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename4 ': v(z)'])
plot(z4, v4, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


figure('NumberTitle', 'Off', 'Name', [filename5 ': z(t)'])
plot(t5, z5, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename5 ': v(t)'])
plot(t5, v5, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename5 ': v(z)'])
plot(z5, v5, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


figure('NumberTitle', 'Off', 'Name', [filename6 ': z(t)'])
plot(t6, z6, '-')
xlabel('t [s]')
ylabel('z [m]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename6 ': v(t)'])
plot(t6, v6, '-')
xlabel('t [s]')
ylabel('v [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', [filename6 ': v(z)'])
plot(z6, v6, '-')
xlabel('z [m]')
ylabel('v [m/s]')
grid on


% Voici un exemple pour les etudes de convergences:
 deltat = [10/200 10/400 10/800 10/1600 10/3200 10/6400];
 zfin = [z1(201) z2(401) z3(801) z4(1601) z5(3201) z6(6401)];
 figure('NumberTitle', 'Off', 'Name', ['Étude de convergence (position, avec atmosphère)'])
 plot(deltat, zfin, '+')
 xlabel('\Deltat [s]')
 ylabel('z(t_{fin}) [m]')
 grid on

 deltat = [10/200 10/400 10/800 10/1600 10/3200 10/6400];
 zfin = [v1(201) v2(401) v3(801) v4(1601) v5(3201) v6(6401)];
 figure('NumberTitle', 'Off', 'Name', ['Étude de convergence (vitesse, avec atmosphère)'])
 plot(deltat, zfin, '+')
 xlabel('\Deltat [s]')
 ylabel('v(t_{fin}) [m/s]')
 grid on
