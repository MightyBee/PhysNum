#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "ConfigFile.tpp"

using namespace std;

// Resolution d'un systeme d'equations lineaires par elimination de Gauss-Jordan:
template <class T>
vector<T> solve(vector<T> const& diag,
                vector<T> const& lower,
                vector<T> const& upper,
                vector<T> const& rhs)
{
  vector<T> solution(diag.size());
  vector<T> new_diag(diag);
  vector<T> new_rhs(rhs);

  for(int i(1); i<diag.size(); ++i)
  {
    double pivot = lower[i-1]/new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  for(int i(diag.size()-2); i>=0; --i)
    solution[i] = (new_rhs[i] - upper[i]*solution[i+1]) / new_diag[i];

  return solution;
}


// Classe pour epsilon_r(r)
class Epsilonr {
public:
  Epsilonr(bool const& trivial_, double const& b_, double const& c_)
    : b(b_), R(c_), trivial(trivial_) {};

  inline double operator()(double const& r, bool const& left) {
  // Le booleen "left" indique s'il faut prendre la limite a gauche ou a droite en cas de discontinuite
    double eps(1e-12*b);
    if(trivial or r<=b-eps or (abs(r-b)<=eps and left))
      return 1.0;
    else
      return 8.0 - 6.0*(r-b)/(R-b);
  }

private:
  double b, R;
  bool trivial;
};

// Classe pour rho_lib(r)/epsilon_0
class Rho_lib {
public:
  Rho_lib(bool const& trivial_, double const& b_, double const& a0_)
    : b(b_), a0(a0_), trivial(trivial_) {};

  inline double operator()(double const& r) {
    if(trivial or r>b)
      return 1.0;
    else
      return a0*(1.0-pow(r/b,2));
  }

private:
  double b, a0;
  bool trivial;
};

int main(int argc, char* argv[])
{

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice6 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice6 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Fichier de sortie :
  string output = configFile.get<string>("output");

  // Domaine :
  const double b(configFile.get<double>("b"));
  const double R(configFile.get<double>("R"));

  // Conditions aux bords :
  const double V0(configFile.get<double>("V0"));

  // Instanciation des objets :
  Epsilonr epsilonr(configFile.get<bool>("trivial"), b, R);
  Rho_lib rho_lib(configFile.get<bool>("trivial"), b, configFile.get<double>("a0"));

  // Discretisation du domaine :
  int N1 = configFile.get<int>("N1");
  int N2 = configFile.get<int>("N2");
  int ninters = N1 + N2;
  int npoints = ninters + 1;
  double h1 = b/N1;
  double h2 = (R-b)/N2;
  vector<double> r(npoints);
  for(int i(0); i<N1; ++i)
    r[i] = i*h1;
  for(int i(0); i<=N2; ++i)
    r[N1+i] = b + i*h2;
  vector<double> h(ninters);
  for(int i(0); i<ninters; ++i)
    h[i] = r[i+1] - r[i];

  vector<double> diag(npoints,0.);  // Diagonale
  vector<double> lower(ninters,0.); // Diagonale inferieure
  vector<double> upper(ninters,0.); // Diagonale superieure
  vector<double> rhs(npoints,0.);   // Membre de droite

  // TODO: Assemblage des elements de la matrice et du membre de droite



  // TODO: Condition au bord:



  // Resolution:
  vector<double> phi(solve(diag,lower,upper,rhs));

  // Export des resultats:
  // 1. phi
  ofstream ofs((output+"_phi.out").c_str());
  ofs.precision(15);
  for(int i(0); i<npoints; ++i)
    ofs << r[i] << " " << phi[i] << endl;
  ofs.close();

  // 2. E_r et D_r
  vector<double> rmid(ninters);
  vector<double> Er(ninters);
  vector<double> Dr(ninters);
  for(int i(0); i<ninters; ++i)
  {
    rmid[i] = 0.5*r[i] + 0.5*r[i+1];
    // TODO: Calculer E_r et D_r/epsilon_0 au milieu des intervalles
    Er[i] = 0.;
    Dr[i] = 0.;
  }
  ofs.open((output+"_Er_Dr.out").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters; ++i)
    ofs << rmid[i] << " " << Er[i] << " " << Dr[i] << endl;
  ofs.close();

  // 3. rho_lib, div(E_r) et div(D_r)
  vector<double> rmidmid(ninters-1);
  vector<double> div_Er(ninters-1);
  vector<double> div_Dr(ninters-1);
  for(int i(0); i<ninters-1; ++i)
  {
    rmidmid[i] = 0.5*rmid[i] + 0.5*rmid[i+1];
    // TODO: Calculer div(E_r) et div(D_r)/epsilon_0 au milieu des milieu des intervalles
    div_Er[i] = 0.;
    div_Dr[i] = 0.;
  }
  ofs.open((output+"_rholib_divEr_divDr.out").c_str());
  ofs.precision(15);
  for(int i(0); i<ninters-1; ++i)
    ofs << rmidmid[i] << " " << rho_lib(rmidmid[i]) << " " << div_Er[i] << " " << div_Dr[i] << endl;
  ofs.close();

  return 0;
}

