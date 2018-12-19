G=6.674e-11;
rho0=0;
tFin=2.5*24*3600;
rL=384748000;
mT=5.972e24;
mL=7.3477e22;
rG=mL*rL/(mT+mL);
rT=0-rG;
rL=rL-rG;
rA=314159000-rG;
omega=sqrt(mL/abs(rT)*G/(rL-rT)^2);

RT=6378100;
RL=1737000;

figure
plot(0,0,'x')
hold on
plot(rT+RT*cos(linspace(0,2*pi,10000)),RT*sin(linspace(0,2*pi,10000)))
hold on
plot(rL+RL*cos(linspace(0,2*pi,10000)),RL*sin(linspace(0,2*pi,10000)))
hold off
grid on
axis equal


