#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

class Engine
{

private:
  double t, tfin; // Temps courant et temps final
  unsigned int nsteps; // Nombre d'iterations
  double B0, Kappa; // Intensite et gradient du champ magnetique
  double E;   // Intensite du champ electrique
  double q;  // Charge de la particule
  double x0, y0, vx0, vy0;  // Position et vitesse initiales de la particle
  unsigned int sampling; // Nombre d'iterations entre chaque ecriture des diagnostics
  unsigned int last; // Nombre d'iterations depuis la derniere ecriture des diagnostics
  ofstream *outputFile; // Pointeur vers le fichier de sortie

  // Ecriture des diagnostics
  void printOut(bool force)
  {
    // Ecriture tous les [sampling] pas de temps, sauf si force est vrai
    if((!force && last>=sampling) || (force && last!=1))
    {
      double energy = 0.5*m*(vx*vx+vy*vy)-q*E*y; // TODO: Completer l'expression de l'energie
      double mu = 0.5*m*(vx*vx+vy*vy)/B(x) ; // TODO: Completer l'expression du moment magnetique
      *outputFile << t << " " << x << " " << y << " " << vx << " " << vy << " " << energy << " " << mu << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // Iteration temporelle, a definir au niveau des classes filles
  virtual void step()=0;

  // Champ magnetique variable
  double B(double const& x) const
  {
    return B0 * (1. + Kappa*x);
  }

protected:

  double ax(double vy, double x){
    return (q*B(x)/m)*vy;
  }
  double ay(double vx, double x){
    return q/m*(E-B(x)*vx);
  }

  double m; // Masse de la particule
  double dt; // Pas de temps
  double x, y, vx, vy;  // Position et vitesse de la particle

public:

  // Constructeur
  Engine(ConfigFile configFile)
  {
    // Stockage des parametres de simulation dans les attributs de la classe
    tfin     = configFile.get<double>("tfin");
    nsteps   = configFile.get<unsigned int>("nsteps");
    dt       = tfin / nsteps;
    m        = configFile.get<double>("m");
    q        = configFile.get<double>("q");
    B0       = configFile.get<double>("B0");
    Kappa    = configFile.get<double>("Kappa");
    E        = configFile.get<double>("E");
    x0       = configFile.get<double>("x0");
    y0       = configFile.get<double>("y0");
    vx0      = configFile.get<double>("vx0");
    vy0      = configFile.get<double>("vy0");
    sampling = configFile.get<unsigned int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
  };

  // Destructeur virtuel
  virtual ~Engine()
  {
    outputFile->close();
    delete outputFile;
  };

  // Simulation complete
  void run()
  {
    t = 0.;
    x = x0;
    y = y0;
    vx = vx0;
    vy = vy0;
    last = 0;
    printOut(true);
    for(unsigned int i(0); i<nsteps; ++i)
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

class EngineEuler: public Engine
{
public:
  EngineEuler(ConfigFile configFile): Engine(configFile) {}

  void step()
  {
    // TODO: Mettre a jour x, y, vx, vy avec le schema d'Euler
    double aX(ax(vy,x)), aY(ay(vx,x));
    vx = vx + aX*dt;
    vy = vy + aY*dt;
    x  = x  + vx*dt;
    y  = y  + vy*dt;
  }
};

class EngineEulerCromer: public Engine
{
public:
  EngineEulerCromer(ConfigFile configFile): Engine(configFile) {}

  void step()
  {
    // TODO: Mettre a jour x, y, vx, vy avec le schema d'Euler-Cromer
    vx = vx + ax(vy,x)*dt;
    vy = vy + ay(vx,x)*dt;
    x  = x  + vx  *dt;
    y  = y  + vy  *dt;
  }
};

class EngineVerlet: public Engine
{
public:
  EngineVerlet(ConfigFile configFile): Engine(configFile) {}

  void step()
  {
    x  = x  + vx*0.5*dt;
    y  = y  + vy*0.5*dt;
    vx = vx + ax(vy,x)*dt/2.0;
    vy = vy + ay(vx,x)*dt/2.0;
    vy = vy + ay(vx,x)*dt/2.0;
    vx = vx + ax(vy,x)*dt/2.0;
    y  = y  + vy*0.5*dt;
    x  = x  + vx*0.5*dt;
  }
};

class EngineRungeKutta2: public Engine
{
public:
  EngineRungeKutta2(ConfigFile configFile): Engine(configFile) {}

  void step()
  {
    // TODO: Mettre a jour x, y, vx, vy avec le schema de Runge-Kutta d'ordre 2
    double kvx1 = dt * ax(vy,x);
    double kvy1 = dt * ay(vx,x);
    double kx1  = dt * vx;
    double ky1  = dt * vy;

    double vxNew = vx + 0.5*kvx1;
    double vyNew = vy + 0.5*kvy1;
    double xNew  = x  + 0.5*kx1;
    double yNew  = y  + 0.5*ky1;

    double kvx2 = dt * ax(vyNew, xNew);
    double kvy2 = dt * ay(vxNew, xNew);
    double kx2  = dt * vxNew;
    double ky2  = dt * vyNew;

    vx = vx + kvx2;
    vy = vy + kvy2;
    x  = x  + kx2 ;
    y  = y  + ky2 ;
  }
};


int main(int argc, char* argv[])
{
  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Schema numerique ("Euler"/"E", "EulerCromer"/"EC" ou "RungeKutta2"/"RK2")
  string schema(configFile.get<string>("schema"));

  Engine* engine;
  if(schema == "Euler" || schema == "E")
  {
    engine = new EngineEuler(configFile);
  }
  else if(schema == "EulerCromer" || schema == "EC")
  {
    engine = new EngineEulerCromer(configFile);
  }
  else if(schema == "Verlet" || schema == "V")
  {
    engine = new EngineVerlet(configFile);
  }
  else if(schema == "RungeKutta2" || schema == "RK2")
  {
    engine = new EngineRungeKutta2(configFile);
  }
  else
  {
    cerr << "Schema inconnu" << endl;
    return -1;
  }

  engine->run();

  delete engine;
  cout << "Fin de la simulation." << endl;
  return 0;
}
