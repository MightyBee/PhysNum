% Nom du fichier d'output a analyser
filename = 'Euler.out';

% Chargement des donnees
output = load(filename);

% Extraction des quantites d'interet
t = output(:,1);
x = output(:,2);
y = output(:,3);
vx = output(:,4);
vy = output(:,5);
energy = output(:,6);
mu = output(:,7);

clear output

% Figures

figure
subplot(2,3,1)
plot(x,y)
axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')

subplot(2,3,2)
plot(vx,vy)
axis equal
grid on
xlabel('v_x [m/s]')
ylabel('v_y [m/s]')

subplot(2,3,3)
plot(t,x,t,y)
grid on
xlabel('t [s]')
ylabel('x,y [m]')
legend('x','y')

subplot(2,3,4)
plot(t,vx,t,vy)
grid on
xlabel('t [s]')
ylabel('v_x,v_y [m/s]')
legend('v_x','v_y')

subplot(2,3,5)
plot(t,energy)
grid on
xlabel('t [s]')
ylabel('E [J]')

subplot(2,3,6)
plot(t,mu)
grid on
xlabel('t [s]')
ylabel('\mu [J/T]')

