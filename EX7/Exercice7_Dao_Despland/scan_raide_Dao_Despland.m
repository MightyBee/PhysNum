% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%
% Il utilise les arguments du programme (voir ConfigFile.h)
% pour remplacer la valeur d'un parametre du fichier d'input
% par la valeur scannee.
%

%% Parametres %%
%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile (NB: enlever le ./ sous Windows)
executable = 'Exercice7'; % Nom de l'executable (NB: ajouter .exe sous Windows)
input = 'configuration_t.in'; % Nom du fichier d'entree de base
dossier='simulations';

nsimul = 4; % Nombre de simulations a faire

xa0=200000;
xb0=370000;
d0=xb0-xa0;

d=round(linspace(d0,0.1*d0, nsimul));

xa=xa0+(d0-d)/2;
xb=xb0-(d0-d)/2;

%% Simulations %%
%%%%%%%%%%%%%%%%%

output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for i = 1:nsimul
    output{i} = sprintf('%s/xa=%d_xb=%d',dossier,xa(i),xb(i));
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s xa=%d xb=%d output=%s', repertoire, executable, input, xa(i), xb(i), output{i});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%

% Parcours des resultats de toutes les simulations

x=cell(nsimul,1);
t=cell(nsimul,1);
f=cell(nsimul,1);

for i = 1:nsimul % Parcours des resultats de toutes les simulations
    data = load([output{i} '_u.out']);
    x{i} = data(:,1);
    u = data(:,2);
    data = load([output{i} '_f.out']);
    t{i} = data(1:end,1);
    f{i} = data(1:end,2:end);
end


%% Figures %%
%%%%%%%%%%%%%

% n=nsimul;
% while 1
%     div=divisors(n)
%     if length(div)==2
%         n=n+1;
%     else
%         n1=div(round(length(div)/2))
%         n2=n/n1
%         if n2/n1<2
%             break;
%         else
%             n=n+1;
%         end
%     end
% end
%     
% 
% 
% 
% pente=pi/2*7980./d;
% fig1=figure;%('Position',[50,50,600,450]);
% for i = 1:nsimul % Parcours des resultats de toutes les simulations
%     subplot(n1,n2,i)
%     pcolor(x{i},t{i},f{i})
%     shading interp
%     colormap jet
%     c = colorbar;
%     xlabel('x [m]')
%     ylabel('t [s]')
%     ylabel(c,'f(x,t) [m]')
%     title(sprintf('Pente_{max} : %d',pente(i)))
% end
% print(fig1,'figures/scan_raide', '-depsc');


fig2=figure('Position',[50,50,600,450]);
pcolor(x{nsimul},t{nsimul},f{nsimul})
shading interp
colormap jet
c = colorbar;
xlabel('$x \ \rm [m]$','Interpreter','Latex')
ylabel('$t \ \rm [s]$','Interpreter','Latex')
ylabel(c,'$f(x,t) \ \rm [m]$','Interpreter','Latex')
title(sprintf('Pente_{max} : %d',pente(nsimul)),'FontSize',8)
grid on, box on
set(gca,'FontSize',25)
pos=get(gca,'position');  % retrieve the current values
pos(3)=0.9*pos(3);        % try reducing width 10%
set(gca,'position',pos);  % write the new values
print(fig2,'figures/zoom_raide', '-depsc');