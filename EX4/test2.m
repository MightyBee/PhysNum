output = load('simulations/h.out');

t = output(:,1);
x1 = output(:,2);
y1 = output(:,3);
x2 = output(:,8);
y2 = output(:,9);
x3 = output(:,14);
y3 = output(:,15);

figure
plot(x1,y1,x2,y2,x3,y3)
grid on