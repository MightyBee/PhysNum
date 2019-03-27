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

double energie(vector<double> const& f, double const& dx)
{
  return 0.;
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
  double L       = configFile.get<double>("L");
  int N          = configFile.get<int>("Npoints");
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

  double dx = L / (N-1);
  double dt = CFL * dx / sqrt(u2->max());
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non
  string schema = configFile.get<string>("schema");
  if(schema != "A" && schema != "B" && schema != "C")
  {
    cerr << "Merci de choisir schema=""A"", ""B"" ou ""C""." << endl;
    return -1;
  }

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

  ofstream fichier_E((output + "_E.out").c_str());
  fichier_E.precision(15);

  ofstream fichier_u((output + "_u.out").c_str());
  fichier_u.precision(15);
  for(double x(0.); x<=L+.5*dx; x+=dx)
    fichier_u << x << " " << sqrt((*u2)(x)) << endl;
  fichier_u.close();


  // Initialisation des tableaux du schema numerique :
  vector<double> fpast(N), fnow(N), fnext(N);

  for(int i(0); i<N; ++i)
  {
    fpast[i] = 0.;
    fnow[i] = 0.;
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
      fichier_E << t << " " << energie(fnow,dx) << endl;
    }
    ++stride;

    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
      if (schema == "A")
	fnext[i] = 0.;// TODO : Completer le schema A
      else if(schema == "B")
        fnext[i] = 0.; // TODO : Completer le schema B
      else if(schema=="C")
        fnext[i] = 0.; // TODO : Completer le schema C
   // Note : La syntaxe pour evaluer u^2 au point x est (*u2)(x)
    }

    // Conditions aux bords :
    switch(cb_gauche)
    {
      case fixe:
        fnext[0] = 0.; // TODO : Completer la condition au bord gauche fixe
        break;

      case libre:
        fnext[0] = 0.; // TODO : Completer la condition au bord gauche libre
        break;

      case harmonique:
        fnext[0] = 0.; // TODO : Completer la condition au bord gauche harmonique
        break;

      case sortie:
        fnext[0] = 0.; // TODO : Completer la condition au bord gauche "sortie de l'onde"
        break;
    }

    switch(cb_droit)
    {
      case fixe:
      fnext[N-1] = 0; // TODO : Completer la condition au bord droit fixe
        break;

      case libre:
        fnext[N-1] = 0.; // TODO : Completer la condition au bord droit libre
        break;

      case harmonique:
        fnext[N-1] = 0.; // TODO : Completer la condition au bord droit harmonique
        break;

      case sortie:
        fnext[N-1] = 0.; // TODO : Completer la condition au bord droit "sortie de l'onde"
        break;
    }

    // Mise a jour :
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_E << t << " " << energie(fnow,dx) << endl;

  fichier_f.close();
  fichier_E.close();

  return 0;
}
