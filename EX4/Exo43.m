%% ConfigFile %%
%%%%%%%%%%%%%%%%


G=6.674e-11;
rho0=1.2;
tFin=2*24*3600;

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             G         % G
             rho0      % rho0
             7238.2    % lambda
             1000        % dt
             1e-5      % precision
             "true"    % adaptatif
             "deuxCorps.out"   % output
             0      ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);

v0=1200;
r0=314159000;
mT=5.972e24;
h=10000;
RT=6378100;
vMax_th=sqrt(v0^2+2*G*mT*(1/(h+RT)-1/r0));
alpha = pi-asin(vMax_th*(h+RT)/(v0*r0));
vx0=v0*cos(alpha);
vy0=v0*sin(alpha);

rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Terre',        'Lune',         'Apollo13'        };
variables = [0               384748000       r0              % x0
             0               0               0               % y0
             0               0               0               % z0
             0               0               vx0             % vx0
             0               0               vy0             % vy0
             0               0               0               % vz0
             mT              7.3477e-22      5809            % m
             RT              1737000         1.95            % R
             0               0               0.3          ]; % Cx

T1=table(variables(:,1),variables(:,2),variables(:,3),'VariableNames',varNames,'RowNames',rowNames);

config(T0,T1);




%% Parametres à varier %%
%%%%%%%%%%%%%%%%%%%%%%%%%

repertoire = './'; % Chemin d'acces au code compile
executable = 'performance'; % Nom de l'executable 

nsimul = 2; % Nombre de simulations à faire

precision = logspace(1,-9,nsimul); % Valeurs du parametre a scanner


    change_config(0,'adaptatif','true');
   paramstr={"precision"};
    param=precision;
   configfileNb=0;


%% Simulations %%
%%%%%%%%%%%%%%%%%




input=sprintf('configuration%d.in', configfileNb);
output = cell(1, nsimul); % Tableau de cellules contenant le nom des fichiers de sortie
for k = 1:nsimul
    parameter = "";
    for i=1:size(paramstr,1)
      parameter=parameter+sprintf('%s=%.15g ', paramstr{i,1}, param(i,k));
    end
    parameter=strip(parameter);
    output{1,k} = "simulations/"+strrep(parameter, ' ', '_')+".out";
    % Execution du programme en lui envoyant la valeur a scanner en argument
    cmd = sprintf('%s%s %s %d %s configuration0.in 1 output=%s', repertoire, executable, input, size(param,1), parameter, output{1,k});
    disp(cmd)
    system(cmd);
end

%% Analyse %%
%%%%%%%%%%%%%
maxAcc=zeros(1,nsimul);
maxPT=zeros(1,nsimul);
nsteps=ones(1,nsimul);

for i = 1:nsimul % Parcours des resultats de toutes les simulations 
    data = load(output{i}); % Chargement du fichier de sortie de la i-ieme simulation
        t = data(:,1);
        xT = data(:,5);
        yT = data(:,5);
        xA = data(:,17);
        yA = data(:,18);
        vx= data(:,20);
        vy= data(:,21);
        nsteps(i)=size(t,1)-1;
        R = T.Terre(8);
        Acc = data(:,2);
        a=Acc(1);
        for l=1:size(Acc,1)
            if Acc(l)>a
                a=Acc(l);
            end
        end
        maxAcc(i)= a; % 
        PT = -data(:,4);
        v =PT(1);
       for l=1:size(PT,1)
            if PT(l)>v
                v=PT(l);
            end
        end
        maxPT(i)= v; 
end


figure 
angle=linspace(0,2*pi);
    plot(RT*cos(angle),RT*sin(angle),'r')
    hold on
    plot(xT(1),yT(1),'+r',xA,yA)
    hold off
    axis equal
    grid on
  
 figure
 loglog(nsteps , maxAcc,'+')
 grid on
 
 
 figure
 loglog(nsteps, maxPT,'+')
 grid on
