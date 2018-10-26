
filenameA = 'v0_97h20.out';
filenameB = 'v0_11km.out';


dataA = load(filenameA);
dataB = load(filenameB);

dtA = dataA(:,2);
v0A = dataA(:,3);
dtB = dataB(:,2);
v0B = dataB(:,3);



figure('NumberTitle', 'Off', 'Name', ['Convergenve : v0(dt)'])
plot(dtA, v0A, '+')
xlabel('\Delta t [s]')
ylabel('v0 [m/s]')
grid on

figure('NumberTitle', 'Off', 'Name', ['Convergenve : v0(dt)'])
plot(dtB, v0B, '+')
xlabel('\Delta t [s]')
ylabel('v0 [m/s]')
grid on
