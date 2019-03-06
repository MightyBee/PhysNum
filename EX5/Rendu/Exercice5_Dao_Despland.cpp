#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include "ConfigFile.tpp"

using namespace std;

double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2);
size_t index(double const& x, double const& h, bool const& low);

int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice5 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice5 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);



  double choix = configFile.get<int>("choix");

  // Geometrie:
  double L  = configFile.get<double>("L");
  double xa = configFile.get<double>("xa");
  double xb = configFile.get<double>("xb");
  double xc = configFile.get<double>("xc");
  double xd = configFile.get<double>("xd");
  double ya = configFile.get<double>("ya");
  double yb = configFile.get<double>("yb");
  double R1  = configFile.get<double>("R1");
  double R2  = configFile.get<double>("R2");
  double R3  = configFile.get<double>("R3");

  // Temperatures:
  double Tc = configFile.get<double>("Tc");
  double Tf = configFile.get<double>("Tf");
  double Tb = configFile.get<double>("Tb");

  // Physical coefficients
  double kappa = configFile.get<double>("kappa");
  double rho = configFile.get<double>("rho");
  double Ccoef = configFile.get<double>("Ccoef");
  double D = kappa/(rho*Ccoef);

  // Duree de la simulation:
  double tfin = configFile.get<double>("tfin");
  double eps = configFile.get<double>("eps"); // Condition d'arret si etat stationnaire

  // Discretisation:
  size_t N = configFile.get<size_t>("N"); // Nombre d'intervalles dans chaque dimension
  double dt = configFile.get<double>("dt");
  double h = L/N;
  double alpha = D * dt / h / h;

  // Fichiers de sortie:
  string output = configFile.get<string>("output");
  ofstream output_T((output+"_T.out").c_str()); // Temperature au temps final
  ofstream output_P((output+"_P.out").c_str()); // Puissance au cours du temps
  output_T.precision(15);
  output_P.precision(15);

  // Tableaux:
  vector<vector<bool> > flag(N+1,vector<bool>(N+1));
  vector<vector<double> > T(N+1,vector<double>(N+1));

  // Initialisation des tableaux
  if(choix==1){
  for(size_t i(0); i<N+1; i++){
    for(size_t j(0); j<N+1; j++){
      if(i==0 || i==N || j==0 || j==N){
        flag[i][j]=true;
        T[i][j]=Tb;
      }else if(i*h>=xa && i*h<=xb && j*h>=ya && j*h<=yb){
        flag[i][j]=true;
        T[i][j]=Tc;
      }else if(i*h>=xc && i*h<=xd && j*h>=ya && j*h<=yb){
        flag[i][j]=true;
        T[i][j]=Tf;
      }else{
        flag[i][j]=false;
        T[i][j]=Tb;
      }
    }
  }
}else if(choix==2){
  for(size_t i(0); i<N+1; i++){
    for(size_t j(0); j<N+1; j++){
      if(i==0 || i==N || j==0 || j==N){
        flag[i][j]=true;
        T[i][j]=Tb;
      }else if(i*h>=xa && i*h<=xb && j*h>=ya+(i*h-xa)*(yb-ya)/(xb-xa) && j*h<=yb){
        flag[i][j]=true;
        T[i][j]=Tc;
      }else if(i*h>=xc && i*h<=xd && j*h>=ya && j*h<=ya+(i*h-xc)*(yb-ya)/(xd-xc)){
        flag[i][j]=true;
        T[i][j]=Tf;
      }else{
        flag[i][j]=false;
        T[i][j]=Tb;
      }
    }
  }
}else if(choix==3){
  xa=(xa+xb)/2;
  xb=(xc+xd)/2;
  ya=(ya+yb)/2;
  for(size_t i(0); i<N+1; i++){
    for(size_t j(0); j<N+1; j++){
      if(i==0 || i==N || j==0 || j==N){
        flag[i][j]=true;
        T[i][j]=Tb;
      }else if(pow(i*h-xa,2)+pow(j*h-ya,2)<=R1*R1){
        flag[i][j]=true;
        T[i][j]=Tc;
      }else if(pow(i*h-xb,2)+pow(j*h-ya,2)<=R2*R2){
        flag[i][j]=true;
        T[i][j]=Tf;
      }else{
        flag[i][j]=false;
        T[i][j]=Tb;
      }
    }
  }
}else{
  xa=(xa+xd)/2;
  ya=(ya+yb)/2;
for(size_t i(0); i<N+1; i++){
  for(size_t j(0); j<N+1; j++){
    if(i==0 || i==N || j==0 || j==N){
      flag[i][j]=true;
      T[i][j]=Tb;
    }else if(pow(i*h-xa,2)+pow(j*h-ya,2)<=R1*R1){
      flag[i][j]=true;
      T[i][j]=Tf;
    }else if(not (i*h<xa && j*h>ya-0.15*R2 && j*h<ya+0.15*R2) && pow(i*h-xa,2)+pow(j*h-ya,2)>=R2*R2 && pow(i*h-xa,2)+pow(j*h-ya,2)<=R3*R3){
      flag[i][j]=true;
      if(j*h==ya && i*h<xa) T[i][j]=Tf+0.5*Tc;
      else T[i][j]=(Tf+Tc)*0.5+(Tc-Tf)*(2*atan((j*h-ya)/((i*h-xa)+sqrt(pow(i*h-xa,2)+pow(j*h-ya,2)))))/(2*M_PI);
    }else{
      flag[i][j]=false;
      T[i][j]=Tb;
    }
  }
}
}

  // Iterations:
  //////////////////////////////////////
  double MaxdT(0.0);
  double Told(0.0);
  double derivee(0.0);
  // TODO: Modifier la condition de sortie de la boucle temporelle pour tester si l'etat stationnaire est atteint.
  size_t iter(0);
  do{
    MaxdT=0.0;
  //for(size_t iter(0); iter*dt<tfin and (MaxdT>eps or iter==0); ++iter){
    // TODO: Schema a 2 niveaux et calcul de max(|dT/dt|)
    for(size_t i(0); i<N+1;i++){
      for(size_t j(0); j<N+1; j++){
        if(flag[i][j]==false){
          Told=T[i][j];
          T[i][j]=T[i][j]+alpha*(T[i-1][j]+T[i+1][j]-4*T[i][j]+T[i][j+1]+T[i][j-1]);
          derivee=abs((T[i][j]-Told)/dt);
          if(derivee>MaxdT){
            MaxdT=derivee;
          }
        }
      }
    }
    // Diagnostiques:
    output_P << iter*dt << " " << puissance(T, kappa, h, xa, xb, ya, yb)
                        << " " << puissance(T, kappa, h, xc, xd, ya, yb)
                        << " " << puissance(T, kappa, h, xa, xd, ya, yb) << endl;
  //}
  iter++;
  }while(iter*dt<tfin and MaxdT>eps);
  output_P.close();

  if(MaxdT<=eps) cout << "Etat stationnaire après " << iter << " iterations." << endl;

  // Ecriture de la temperature finale:
  for(size_t i(0);i<N+1;++i)
    for(size_t j(0);j<N+1;++j)
      output_T << i*h << " " << j*h << " " << T[i][j] << endl;
  output_T.close();
  return 0;
}

// TODO: Calculer la puissance calorifique emise/recue par le rectangle allant de (x1,y1) a (x2,y2)
double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2)
{
  double P(0.0);

  size_t indexFin(index(x2,h,false));
  size_t ind1(index(y1,h,true));
  size_t ind2(index(y2,h,false));
  for(size_t i(index(x1,h,true));i<=indexFin;i++){
    P+=T[i][ind2+1]-T[i][ind2];
    P-=T[i][ind1]  -T[i][ind1-1];
  }

  indexFin=index(y2,h,false);
  ind1=index(x1,h,true);
  ind2=index(x2,h,false);
  for(size_t j(index(y1,h,true));j<=indexFin;j++){
    P+=T[ind2+1][j]-T[ind2][j];
    P-=T[ind1][j]  -T[ind1-1][j];
  }
  P*=-kappa;
  return P;
}

size_t index(double const& x, double const& h, bool const& low){
  if(low) return (size_t)(x/h-0.5);
  else return (size_t)(x/h+0.5)+1;
}
