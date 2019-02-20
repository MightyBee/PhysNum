% Ce script Matlab automatise la production de resultats
% lorsqu'on doit faire une serie de simulations en
% variant un des parametres d'entree.
%

%% ConfigFile %%
%%%%%%%%%%%%%%%%

tFin=10*365*24*3600;
output="solaire.out";

rowNames  = {'nbCorps','tFin','G','rho0','lambda','dt','precision','adaptatif','output', 'sampling'};
varNames  = {'classique'}; % nom
variables = [3         % nbCorps
             tFin      % tFin
             6.674e-11 % G
             0         % rho0
             7238.2    % lambda
             50        % dt
             1     % precision
             "true"   % adaptatif
             output    % output
             1     ]; % sampling

T0=table(variables,'VariableNames',varNames,'RowNames',rowNames);


rowNames  = {'x0','y0','z0','vx0','vy0','vz0','m','R','Cx'};
varNames  = {'Soleil',  'Mercure', 'Venus',   'Terre',   'Mars',    'Jupiter',  'Saturne',  'Uranus',    'Neptune'  };
variables = [[0          46001272   107476259  147098074  206644545  740520000   1349823615  2734998229   4452940833]*1e3 % x0
              0          0          0          0          0          0           0           0           0          % y0
              0          0          0          0          0          0           0           0           0          % z0
              0          0          0          0          0          0           0           0           0          % vx0
              0          58980      35260      30287      26499      13720       10183       7128        5479       % vy0
              0          0          0          0          0          0           0           0           0          % vz0
              1.9891e30  3.3011e23  4.8685e24  5.9736e24  641.85e21  1.8986e27   568.46e24   8.6810e25   102.43e24  % m
              0          0          0          0          0          0           0           0           0          % R
              0          0          0          0          0          0           0           0           0        ];% Cx

          
T1=table(variables(:,1),variables(:,2),variables(:,3),variables(:,4),variables(:,5),variables(:,6),variables(:,7),variables(:,8),variables(:,9),'VariableNames',varNames,'RowNames',rowNames);
nbCorps=size(variables,2);

config(T0,T1);
change_config(0,'nbCorps',nbCorps);

%% simulation %%

cmd = "./Exercice4";
disp(cmd)
system(cmd);


%% Analyse %%
%%%%%%%%%%%%%


x=cell(1,nbCorps);
y=cell(1,nbCorps);

% 1 2   3  4     5  6  7  8   9   10     11 12 13  14  15  16     17 18 19  20  21  22
% t acc en Pt    x1 y1 z1 vx1 vy1 vz1    x2 y2 z2  vx2 vy2 vz2  


data = load("simulations/"+output); % Chargement du fichier de sortie de la i-ieme simulation
t = data(:,1);
for i=1:nbCorps
   x{i}=data(:,5+(i-1)*6); 
   y{i}=data(:,6+(i-1)*6);
end

%% Figures %%

figure 
plot(x{1},y{1},'p')
hold on
for i=2:nbCorps
   plot(x{i},y{i})
   hold on
end
hold off
axis equal
grid on