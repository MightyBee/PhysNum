#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.tpp"

using namespace std;
typedef vector<complex<double> > vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// ou A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution)
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // forward elimination
  for(int i(1); i<diag.size(); ++i)
  {
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution.resize(diag.size());

  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  // backward substitution
  for(int i = diag.size() - 2; i >= 0; --i)
  {
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) :
double V(double const& x, double const& omega, double const& delta)
{
  return .5*omega*omega*min((x-delta)*(x-delta),(x+delta)*(x+delta));
}

double trapezes(vec_cmplx const& f, double const& dx);

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob calcule la probabilite de trouver la particule entre les points nL.dx et nR.dx,
//  - E calcule son energie,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carre moyenne,
//  - pmoy calcule sa quantite de mouvement moyenne,
//  - p2moy calcule sa quantite de mouvement au carre moyenne.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.


int main(int argc,char **argv)
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres physiques :
  double hbar    = 1.;
  double m       = 1.;
  double tfin    = configFile.get<double>("tfin");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  double omega   = configFile.get<double>("omega");
  double delta   = configFile.get<double>("delta");
  double x0      = configFile.get<double>("x0");
  double k0      = 2. * M_PI * configFile.get<int>("n") / (xR-xL);
  double sigma0  = configFile.get<double>("sigma_norm") * (xR-xL);

  // Parametres numeriques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // TODO: initialiser le paquet d'onde, equation (4.109) du cours
  for(int i(0); i<Npoints; ++i)
    psi[i] = exp(complex<double>(0.0,1.0)*k0*x[i])*exp(-(x[i]-x0)*(x[i]-x0)/(2*sigma_norm*sigma_norm));
  // Modifications des valeurs aux bords :
  psi[0] = complex<double> (0.,0.);
  psi[Npoints-1] = complex<double> (0.,0.);
  // Normalisation :
  psi = normalize(psi, dx);

  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'equation (4.90)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'equation (4.90)

  complex<double> a, b;

  // TODO: calculer les elements des matrices A, B et H.
  // Ces matrices sont stockees sous forme tridiagonale, d:diagonale, c et a: diagonales superieures et inferieures
  for(size_t i(0); i<Npoints; i++){
    dH[i]=1.0/(dx*dx)+V(i*dx,omega,delta);
    dA[i]=complex<double>(1.0,0.0)+complex<double>(0.0,1.0)*dt*0.5*dH[i];
    dB[i]=complex<double>(1.0,0.0)-complex<double>(0.0,1.0)*dt*0.5*dH[i];
  }

  for(size_t i(0); i<Ninters; i++){
    aH[i]=-1.0/(2.0*dx*dx);
    cH[i]=-1.0/(2.0*dx*dx);
    aA[i]=complex<double>(0.0,1.0)*dt*0.5*aH[i];
    cA[i]=complex<double>(0.0,1.0)*dt*0.5*cH[i];
    aB[i]=-complex<double>(0.0,1.0)*dt*0.5*aH[i];
    cB[i]=-complex<double>(0.0,1.0)*dt*0.5*cH[i];
  }

  // Conditions aux limites: psi nulle aux deux bords
  // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
  dA[0]=1;
  cA[0]=0;
  dB[0]=1;
  cB[0]=0;
  dA[Npoints-1]=1;
  aA[Ninters-1]=0;
  dB[Npoints-1]=1;
  aB[Ninters-1]=0;

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_potentiel((output + "_pot.out").c_str());
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega, delta) << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((output + "_psi2.out").c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((output + "_obs.out").c_str());
  fichier_observables.precision(15);

  // Boucle temporelle :
  double t;
  for(t=0.; t+dt/2.<tfin; t+=dt)
  {
    // Ecriture de |psi|^2 :
    for(int i(0); i<Npoints; ++i)
      fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";
    fichier_psi << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,Ninters*xL/(xL-xR),dx) << " "       // probabilite que la particule soit en x < 0
                        << prob(psi,Ninters*xL/(xL-xR),Ninters,dx) << " " // probabilite que la particule soit en x > 0
                        << E(psi,dH,aH,cH,dx) << " "                      // Energie
                        << xmoy(psi,x,dx) << " "                          // Position moyenne
                        << x2moy(psi,x,dx) << " "                         // Position^2 moyenne
                        << pmoy(psi,dx) << " "                            // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << endl;                         // (Quantite de mouvement)^2 moyenne

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints,0.);

    // Multiplication psi_tmp = B * psi :
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i)
    {
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }

    // Resolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);

  } // Fin de la boucle temporelle

  for(int i(0); i<Npoints; ++i)
    fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";

  fichier_observables << t << " "
                      << prob(psi,0,Ninters*xL/(xL-xR),dx) << " "
                      << prob(psi,Ninters*xL/(xL-xR),Ninters,dx) << " "
                      << E(psi,dH,aH,cH,dx) << " "
                      << xmoy(psi,x,dx) << " "
                      << x2moy(psi,x,dx) << " "
                      << pmoy(psi,dx) << " "
                      << p2moy(psi,dx) << endl;

  fichier_observables.close();
  fichier_psi.close();

}



double trapezes(vec_cmplx const& f, double const& dx){
  double retour(0.0);
  for(size_t i(0); i<f.size()-1;i++){
    retour+=abs((f[i]+f[i+1])*0.5*dx);
  }
  return retour;
}

double prob(vec_cmplx const& psi, int nL, int nR, double dx)
{
  // TODO: calculer la probabilite de trouver la particule entre les points nL.dx et nR.dx
  size_t n(abs(nR-nL)+1);
  vec_cmplx psi_rogne(n);
  for(size_t i(0); i<n; i++){
    psi_rogne[i]=psi[i+nL];
  }
  return trapezes(psi_rogne,dx);
}


double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx)
{
  vec_cmplx psi_tmp(psi.size());

  // TODO: calculer la moyenne de l'Hamiltonien

  // H(psi)
  // On utilise la matrice H calculée plus haut
  psi_tmp[0]=diagH[0]*psi[0]+upperH[0]*psi[1];
  for(size_t i(1); i<psi.size()-1; i++){
    psi_tmp[i]=lowerH[i-1]*psi[i-1]+diagH[i]*psi[i]+upperH[i]*psi[i+1];
  }
  psi_tmp[psi.size()-1]=lowerH[psi.size()-2]*psi[psi.size()-2]+diagH[psi.size()-1]*psi[psi.size()-1];
  // Integrale de psi* H(psi) dx
  for(size_t i(0); i<psi.size(); i++){
    psi_tmp[i]=conj(psi[i])*psi_tmp[i];
  }
  return trapezes(psi_tmp,dx);
}


double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  vec_cmplx psi_tmp(psi.size());
  for(size_t i(0); i<psi.size(); i++){
    psi_tmp[i]=conj(psi[i])*x[i]*psi[i];
  }
  return trapezes(psi_tmp,dx);
}


double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  vec_cmplx psi_tmp(psi.size());
  for(size_t i(0); i<psi.size(); i++){
    psi_tmp[i]=conj(psi[i])*x[i]*x[i]*psi[i];
  }
  return trapezes(psi_tmp,dx);
}


double pmoy(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_tmp(psi.size());
  psi_tmp[0]=complex<double>(0.0,-1.0)*conj(psi[0])*(psi[1]-psi[0])/(dx);
  for(size_t i(1); i<psi.size()-1; i++){
    psi_tmp[i]=complex<double>(0.0,-1.0)*conj(psi[i])*(psi[i+1]-psi[i-1])/(2*dx);
  }
  psi_tmp[psi.size()-1]=complex<double>(0.0,-1.0)*conj(psi[psi.size()-1])*(psi[psi.size()-1]-psi[psi.size()-2])/(dx);
  return trapezes(psi_tmp,dx);
}


double p2moy(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_tmp(psi.size());
  psi_tmp[0]=complex<double>(0.0,0.0);
  for(size_t i(1); i<psi.size()-1; i++){
    psi_tmp[i]=-conj(psi[i])*(psi[i+1]-2.0*psi[i]+psi[i-1])/(dx*dx);
  }
  psi_tmp[psi.size()-1]=complex<double>(0.0,0.0);
  return trapezes(psi_tmp,dx);
}

vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi,0,psi.size()-1,dx));
  for(unsigned int i(0); i<psi.size(); ++i)
    psi_norm[i] = psi[i]/norm;
  return psi_norm;
}
