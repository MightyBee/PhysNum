%function Analyse(filename)
filename='output';
data = load([ filename '_Er_Dr.out']);
rmid = data(:,1);
Er = data(:,2);
Dr = data(:,3);
data = load([ filename '_phi.out']);
r = data(:,1);
phi = data(:,2);
data = load([ filename '_rholib_divEr_divDr.out']);
rmidmid = data(:,1);
rholib = data(:,2);
divEr = data(:,3);
divDr = data(:,4);

fig1=figure('Position',[50,50,600,450]);
plot(r,phi)
xlabel('r [m]')
ylabel('\phi [V]')
grid on, box on
set(gca,'FontSize',20)
if(strcmp(filename,'trivial'))
    hold on
    plot(r,(-r.^2.+r(end)^2)/4,'--')
    lgd=legend('Valeurs numériques', 'Solution analytique');
    set(lgd,'fontsize',14,'Location','southwest');
end
print(fig1,['figures/phi_' filename], '-depsc');

fig2=figure('Position',[50,50,600,450]);
hold on
plot(rmid,Er)
xlabel('r [m]')
ylabel('E_r [V/m]')
grid on, box on
set(gca,'FontSize',20)
if(strcmp(filename,'trivial'))
    hold on
    plot(rmid,rmid/2,'--')
    lgd=legend('Valeurs numériques', 'Solution analytique');
    set(lgd,'fontsize',14,'Location','northwest');
end
print(fig2,['figures/Er_' filename], '-depsc');

fig3=figure('Position',[50,50,600,450]);
hold on
plot(rmid,Dr)
xlabel('r [m]')
ylabel('D_r/\epsilon_0 [V/m]')
grid on, box on
set(gca,'FontSize',20)
print(fig3,['figures/Dr_' filename], '-depsc');

fig4=figure('Position',[50,50,600,450]);
hold on
plot(rmidmid,rholib,'DisplayName','\rho_{lib}/\epsilon_0')
plot(rmidmid,divDr,'--','DisplayName','div(D_r)/\epsilon_0')
plot(rmidmid,divEr-divDr,'DisplayName','\rho_{pol}/\epsilon_0')
xlabel('r')
ylabel('\rho/\epsilon_0 [V/m^2]')
grid on, box on
set(gca,'FontSize',20)
lgd=legend('show');
set(lgd,'fontsize',14,'Location','southeast');
print(fig4,['figures/rho_' filename], '-depsc');


