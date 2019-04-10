#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"

using namespace std;

//
// Objets pour la fonction u^2(x)
//
class U2 {
public:
  // Methodes virtuelles pures => classe abstraite
  virtual double operator()(double const& x) const = 0; // Evalue u^2 au point x
  virtual double max() const = 0; // Renvoie max(u^2(x))
};

class U2_const: public U2 {
public:
  // Pas de constructeur par defaut => on force a specifier une valeur
  U2_const(ConfigFile const& configFile) :
  U2(), u2(pow(configFile.get<double>("u"),2))
  {}

  // Definition des methodes virtuelles pures :
  double operator()(double const& x) const
  {
    return u2;
  }

  double max() const
  {
    return u2;
  }

private:
  double u2;
};

class U2_tsunami: public U2 {
public:

  U2_tsunami(ConfigFile const& configFile) :
  U2(),
  g(configFile.get<double>("g")),
  h_ocean(configFile.get<double>("h_ocean")),
  h_recif(configFile.get<double>("h_recif")),
  xa(configFile.get<double>("xa")),
  xb(configFile.get<double>("xb")),
  xc(configFile.get<double>("xc")),
  xd(configFile.get<double>("xd"))
  {}

  double operator()(double const& x) const
  {
    if(x < xa)
      return g * h_ocean;
    else if(x < xb)
      return g * (h_ocean + (h_recif-h_ocean)*pow(sin(.5*M_PI*(x-xa)/(xb-xa)),2));
    else if(x < xc)
      return g * h_recif;
    else if(x < xd)
      return g * (h_recif - (h_recif-h_ocean)*pow(sin(.5*M_PI*(xc-x)/(xc-xd)),2));
    else
      return g * h_ocean;
  }

  double max() const
  {
    return g * std::max(h_ocean, h_recif);
  }

private:
  double g, h_ocean, h_recif, xa, xb, xc, xd;
};

//
// TODO : Calcul de l'energie de l'onde
//

double energie(vector<double> const& f, double const& dh)
{
  double E(0.0);
  for(size_t i(0); i<f.size()-1; i++){
    E+=(f[i]*f[i]+f[i+1]*f[i+1])*0.5*dh;
  }
  return E;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());

  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";

  if(len > 0)
    o << v[len-1];

  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  double Lx      = configFile.get<double>("Lx");
  double Ly    = configFile.get<double>("Ly");
  int Nx         = configFile.get<int>("Npoints");
  double CFL     = configFile.get<double>("CFL");
  string type_u2 = configFile.get<string>("type_u2");

  U2* u2;
  if(type_u2 == "const")
    u2 = new U2_const(configFile);
  else if(type_u2 == "tsunami")
    u2 = new U2_tsunami(configFile);
  else
  {
    cerr << "Merci de choisir type_u2=""const"" ou ""tsunami""." << endl;
    return -1;
  }


  double dh(Lx/(Nx-1));
  int Ny((int)ceil(Ly/dh)+1);
  double dt = CFL * dh / sqrt(u2->max());
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non

  // Conditions aux bords (les strings sont converties en valeurs numeriques a l'aide d'un enumerateur) :
  typedef enum{fixe,libre,harmonique,sortie} Cond_bord;
  Cond_bord cb_gauche, cb_droit;

  string cb = configFile.get<string>("cb_gauche");
  if(cb == "fixe")
    cb_gauche = fixe;
  else if(cb == "libre")
    cb_gauche = libre;
  else if(cb == "harmonique")
    cb_gauche = harmonique;
  else if(cb == "sortie")
    cb_gauche = sortie;
  else
  {
    cerr << "Merci de choisir cb_gauche=""fixe"", ""libre"", ""harmonique"", ou ""sortie""." << endl;
    return -1;
  }

  cb = configFile.get<string>("cb_droit");
  if(cb == "fixe")
    cb_droit = fixe;
  else if(cb == "libre")
    cb_droit = libre;
  else if(cb == "harmonique")
    cb_droit = harmonique;
  else if(cb == "sortie")
    cb_droit = sortie;
  else
  {
    cerr << "Merci de choisir cb_droit=""fixe"", ""libre"", ""harmonique"", ou ""sortie""." << endl;
    return -1;
  }

  double A, omega; // Parametres d'excitation
  if(cb_gauche == harmonique || cb_droit == harmonique)
  {
    A = configFile.get<double>("A");
    omega = configFile.get<double>("omega");
  }


  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_f((output + "_f.out").c_str());
  fichier_f.precision(15);

  ofstream fichier_u((output + "_u.out").c_str());
  fichier_u.precision(15);
  for(double x(0.); x<=Lx+.5*dh; x+=dh)
    fichier_u << x << " " << sqrt((*u2)(x)) << endl;
  fichier_u.close();


  // Initialisation des tableaux du schema numerique :
  vector<vector<double>> fpast(Nx,vector<double>(Ny)), fnow(Nx,vector<double>(Ny)), fnext(Nx,vector<double>(Ny));

  for(int i(0); i<Nx; ++i){
    for(int j(0); j<Ny; ++j){
        fpast[i][j] = 3*exp(-(j*dh-5)*(j*dh-5)-(i*dh-2)*(i*dh-2))+3*exp(-0.8*(j*dh-8)*(j*dh-8)-0.07*(i*dh-5)*(i*dh-5));  cos(sqrt((j*dh-5)*(j*dh-5)+(i*dh-5)*(i*dh-5))*2);
        fnow[i][j]  = 3*exp(-(j*dh-5)*(j*dh-5)-(i*dh-2)*(i*dh-2))+3*exp(-0.8*(j*dh-8)*(j*dh-8)-0.07*(i*dh-5)*(i*dh-5));  cos(sqrt((j*dh-5)*(j*dh-5)+(i*dh-5)*(i*dh-5))*2);
    }
  }


  // Boucle temporelle :
  double t;
  int stride(0);
  int n_stride(configFile.get<int>("n_stride"));
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
    }
    ++stride;

    // Evolution :
    for(int i(1); i<Nx-1; ++i){
      for(int j(1); j<Ny-1; ++j){
        //if((i!=Nx/2 && i!=Nx/2-1 && i!=Nx/2+1) || (j>2.0*Ny/5.0 && j<3.0*Ny/5.0)){
        fnext[i][j] = 2*(1-2*(*u2)(i*dh)*dt*dt/(dh*dh))*fnow[i][j]-fpast[i][j]+(*u2)(i*dh)*dt*dt/(dh*dh)*(fnow[i+1][j]+fnow[i-1][j]+fnow[i][j+1]+fnow[i][j-1]);// TODO : Completer le schema A
      }//}
    }

    // Conditions aux bords :
    switch(cb_gauche)
    {
      case fixe:
        fnext[0] = fnow[0]; // TODO : Completer la condition au bord gauche fixe
        break;

      case libre:
        fnext[0] = fnext[1]; // TODO : Completer la condition au bord gauche libre
        break;

      case harmonique:
        fnext[0] = vector<double>(Ny,A*sin(omega*t));
        /*fnext[Nx-1] = vector<double>(Ny,A*sin(omega*t));
        for(auto& el : fnext){
          el.front()=A*sin(omega*t);
          el.back()=A*sin(omega*t);
        }*/// TODO : Completer la condition au bord gauche harmonique
        break;

      /*case sortie:
        fnext[0] = fnow[0]+sqrt((*u2)(0.5*dh))*dt/dh*(fnow[1]-fnow[0]); // TODO : Completer la condition au bord gauche "sortie de l'onde"
        break;*/
    }

    switch(cb_droit)
    {
      case fixe:
        fnext[Nx-1] = fnow[Nx-1]; // TODO : Completer la condition au bord droit fixe
        break;

      case libre:
        fnext[Nx-1] = fnext[Nx-2]; // TODO : Completer la condition au bord droit libre
        break;

      /*case harmonique:
        fnext[Nx-1] = A*sin(omega*t); // TODO : Completer la condition au bord droit harmonique
        break;

      /*case sortie:
        fnext[Nx-1] = fnow[Nx-1]-sqrt((*u2)((Nx-0.5)*dh))*dt/dh*(fnow[Nx-1]-fnow[Nx-2]);; // TODO : Completer la condition au bord droit "sortie de l'onde"
        break;*/
    }
    /*for(auto& el : fnext){
      el.front()=el[1];
      el.back()=el[Ny-2];
    }*/
    for(size_t i(0); i<Nx; i++){
      fnext[i].front()=fnow[i].front();
      fnext[i].back()=fnow[i].back();
    }
    /*for(size_t j(0); j<Ny; j++){
      if(j<2.0*Ny/5.0 || j>3.0*Ny/5.0){
        fnext[Nx/2][j]=fnow[Nx/2][j];
        fnext[Nx/2-1][j]=fnext[Nx/2-2][j];
        fnext[Nx/2+1][j]=fnext[Nx/2+2][j];
      }
    }*/
    // Mise a jour :
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;

  fichier_f.close();

  return 0;
}
