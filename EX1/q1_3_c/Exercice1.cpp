#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

class Engine
{

private:
  double t, tfin, dt; // Temps courant, temps final et pas de temps
  unsigned int nsteps; // Nombre d'iterations
  double G; // Constante gravitationnelle
  double mT, mL, zT, zL; // Masses et positions de la Terre et de la Lune
  double m, R; // Masse et rayon du projectile
  double z0, v0, v0prime; // Position et vitesse initiales du projectile
  double z, v; // Position et vitesse du projectile au cours du temps
  double rho0, lambda; // Parametres de la densite de l'air
  double Cx; // Coefficient de trainee aerodynamique
  unsigned int sampling; // Nombre d'iterations entre chaque ecriture des diagnostics
  unsigned int last; // Nombre d'iterations depuis la derniere ecriture des diagnostics
  ofstream *outputFile; // Pointeur vers le fichier de sortie

  // Ecriture des diagnostics
  void printOut(bool force)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si force est vrai
    if((!force && last>=sampling) || (force && last!=1))
    {
      *outputFile << nsteps << " " << dt << " " << v0 << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle
  void step()
  {
    double a(Fres()/m);
    v=v+dt*a;
    z=z+dt*v;
  }


public:

  // Constructeur
  Engine(int argc, char* argv[]) {

    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice1 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice1 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");
    nsteps   = configFile.get<unsigned int>("nsteps");
    dt       = tfin / nsteps;
    G        = configFile.get<double>("G");
    mT       = configFile.get<double>("mT");
    mL       = configFile.get<double>("mL");
    zT       = configFile.get<double>("zT");
    zL       = configFile.get<double>("zL");
    m        = configFile.get<double>("m");
    R        = configFile.get<double>("R");
    z0       = configFile.get<double>("z0");
    v0       = configFile.get<double>("v0");
    rho0     = configFile.get<double>("rho0");
    lambda   = configFile.get<double>("lambda");
    Cx       = configFile.get<double>("Cx");
    sampling = configFile.get<unsigned int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur
  ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };


  // algorithme pour trouver la bonne vitesse initiale
  long double findV0(double e=1){
    double v0moins(v0*0.5), v0plus(v0*1.5);
    double vObjectif(11000);
    double zSurfaceLune(zL-1737000);
    for(unsigned int k(0); k<7; k++){
      if(k>0){
        nsteps*=2;
        dt=tfin/nsteps;
        v0moins=v0*0.5;
        v0plus=v0*1.5;
      }
      if(v0plus < 11070){
        v0plus=20000;
      }
      if(v0moins<11070){
        v0moins=11070;
      }
      v0=0.5*(v0moins+v0plus);


      int j(0);
      do{
        j++;
        t = 0;
        z = z0;
        v = v0;
        for(unsigned int i(0); i<nsteps; ++i){
          step();
          t += dt;
        }
        if(z>zSurfaceLune+e){
          v0plus=v0;
          v0=0.5*(v0moins+v0plus);
        }
        if(z<zSurfaceLune-e){
          v0moins=v0;
          v0=0.5*(v0moins+v0plus);
        }
        cout << t << "v0moins : " << v0moins << endl;
        cout << "v0      : " << v0 << endl;
        cout << "v0plus  : " << v0plus << endl << endl;
      }while(j<100 );//and (z<zSurfaceLune-e or z>zSurfaceLune+e));

      /*do{
        t = 0;
        z = z0;
        v = v0;
        for(unsigned int i(0); i<nsteps; ++i){
          step();
          t += dt;
        }
        if(v>vObjectif+e){
          v0plus=v0;
          v0=0.5*(v0moins+v0plus);
        }else if(v<vObjectif-e){
          v0moins=v0;
          v0=0.5*(v0moins+v0plus);
        }
        cout << "v0moins : " << v0moins << endl;
        cout << "v0      : " << v0 << endl;
        cout << "v0plus  : " << v0plus << endl << endl;
      }while(v<vObjectif-e or v>vObjectif+e);*/


      printOut(false);
      cout << nsteps << " steps : v0 = " << v0 << endl;
    }
      return v0;
  }

  double rho(){
    return rho0*exp((z0-z)/lambda);
  }

  double Ft(){
    return -rho()*M_PI*R*R*Cx*v*abs(v)*0.5;
  }

  double FgT(){
    return -m*mT*G/((z-zT)*abs(z-zT));
  }

  double FgL(){
    return -m*mL*G/((z-zL)*abs(z-zL));
  }

  double Fres(){
    return Ft()+FgT()+FgL();
  }

};


int main(int argc, char* argv[])
{
  Engine engine(argc, argv);
  //engine.run();
  cout << fixed << engine.findV0() << endl;
  cout << "Fin de la simulation." << endl;
  return 0;
}
