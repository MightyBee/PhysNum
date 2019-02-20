#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <fstream>
#include "ConfigFile.tpp"

using namespace std;

double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2);

int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice5 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice5 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Geometrie:
  double L  = configFile.get<double>("L");
  double xa = configFile.get<double>("xa");
  double xb = configFile.get<double>("xb");
  double xc = configFile.get<double>("xc");
  double xd = configFile.get<double>("xd");
  double ya = configFile.get<double>("ya");
  double yb = configFile.get<double>("yb");

  // Temperatures:
  double Tc = configFile.get<double>("Tc");
  double Tf = configFile.get<double>("Tf");
  double Tb = configFile.get<double>("Tb");
  double kappa = configFile.get<double>("kappa");

  // Duree de la simulation:
  double tfin = configFile.get<double>("tfin");
  double eps = configFile.get<double>("eps"); // Condition d'arret si etat stationnaire

  // Discretisation:
  int N = configFile.get<int>("N"); // Nombre d'intervalles dans chaque dimension
  double dt = configFile.get<double>("dt");
  double h = L/N;
  double alpha = kappa * dt / h / h;

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

  // Iterations:
  //////////////////////////////////////
  // TODO: Modifier la condition de sortie de la boucle temporelle pour tester si l'etat stationnaire est atteint.
  for(int iter=0; iter*dt<tfin; ++iter)
  {
    // TODO: Schema a 2 niveaux et calcul de max(|dT/dt|)


    // Diagnostiques:
    output_P << iter*dt << " " << puissance(T, kappa, h, xa, xb, ya, yb)
                        << " " << puissance(T, kappa, h, xc, xd, ya, yb)
                        << " " << puissance(T, kappa, h, xa, xd, ya, yb) << endl;
  }
  output_P.close();

  // Ecriture de la temperature finale:
  for(int i(0);i<N+1;++i)
    for(int j(0);j<N+1;++j)
      output_T << i*h << " " << j*h << " " << T[i][j] << endl;
  output_T.close();
  return 0;
}

// TODO: Calculer la puissance calorifique emise/recue par le rectangle allant de (x1,y1) a (x2,y2)
double puissance(vector<vector<double> > const& T, double const& kappa, double const& h, double const& x1, double const& x2, double const& y1, double const& y2)
{
  return 0.;
}
