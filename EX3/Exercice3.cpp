#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double m, g, L;
  double d, Omega, kappa;
  double theta, thetadot;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = 1/2*m*pow(L*thetadot,2)+m*L*g*(1-cos(theta)); // TODO: Evaluer l'energie mecanique
      double pnc = -kappa*pow(L*thetadot,2)-m*pow(Omega,3)*pow(d,2)*sin(Omega*t)*cos(Omega*t); // TODO: Evaluer la puissance des forces non conservatives
      //double pnc = -(kappa*L*thetadot + m*pow(Omega,2)*d*sin(Omega*t)*sin(theta))*(L*thetadot);

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  double a(double x, double v, double t){
    return -(g+d*Omega*Omega/m*sin(Omega*t))*sin(x)/L -kappa/m*v;
  }

  void step()
  {
    double thetaOld(theta);
    theta= theta + thetadot*dt + a(theta,thetadot,t)*dt*dt/2;
    double vTmp(thetadot + a(theta,thetadot,t)*dt/2);
    thetadot=thetadot + (a(thetaOld,vTmp,t)+a(theta,vTmp,t+dt))*dt/2;
  }


public:

  Exercice3(int argc, char* argv[])
  {
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin     = configFile.get<double>("tFin");
    dt       = configFile.get<double>("dt");
    d        = configFile.get<double>("d");
    Omega    = configFile.get<double>("Omega");
    kappa    = configFile.get<double>("kappa");
    m        = configFile.get<double>("m");
    g        = configFile.get<double>("g");
    L        = configFile.get<double>("L");
    theta    = configFile.get<double>("theta0");
    thetadot = configFile.get<double>("thetadot0");
    sampling = configFile.get<int>("sampling");

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
  };

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.;
    last = 0;
    printOut(true);
    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};


int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();
  return 0;
}
