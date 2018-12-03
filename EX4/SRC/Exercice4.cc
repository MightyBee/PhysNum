#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <array>
#include <string>
#include <initializer_list>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template
#include "Vecteur.h"
#include "Systeme.h"
#include "Erreur.h"
using namespace std;


class Exercice4
{

private:
  int nbCorps;
  double t, dt, precision, tFin;
  bool adaptatif;
  Systeme systeme;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << t << " " << systeme << endl;
      last = 1;
    }
    else{
      last++;
    }
  }



public:

  Exercice4(int argc, char* argv[])
  {
    string inputPath0("configuration0.in");
    vector<string> inputPath(0); // Fichier d'input par defaut
    vector<vector<string>> param(0);
    int n(0);
    int k(1);
    while(k<=argc-2){
    	inputPath.push_back(argv[k++]);
      param.push_back(vector<string>(0));
    	n=k+argv[k][0]-'0';
    	while(++k<=n and k<=argc){
    		 param.back().push_back(argv[k]);
    	}
    }

    ConfigFile configFile("config/"+inputPath0);

    for(size_t k(0);k<inputPath.size();k++){
      if(inputPath[k]==inputPath0){
        for(const auto& str : param[k]){
          configFile.process(str);
        }
      }
    }

    // Parametres généraux
    tFin      = configFile.get<double>("tFin");
    dt        = configFile.get<double>("dt");
    precision = configFile.get<double>("precision");
    adaptatif = configFile.get<bool>("adaptatif");
    systeme   = Systeme(configFile.get<double>("G"),configFile.get<double>("rho0"),configFile.get<double>("lambda"));
    nbCorps   = configFile.get<int>("nbCorps");
    sampling  = configFile.get<int>("sampling");
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    // Ajout des corps
    string nom;
    Vecteur P,V;
    double m,R,Cx;
    Corps corps;
    for(int i(1); i<=nbCorps; i++){
      inputPath0="configuration"+to_string(i)+".in"; // Fichier d'input par defaut
      ConfigFile configFile("config/"+inputPath0); // Les parametres sont lus et stockes dans une "map" de strings.
      for(size_t k(0);k<inputPath.size();k++){
        if(inputPath[k]==inputPath0){
          for(const auto& str : param[k]){
            configFile.process(str);
          }
        }
      }
      nom = configFile.get<string>("nom");
      P = Vecteur(configFile.get<double>("x0"),configFile.get<double>("y0"),configFile.get<double>("z0"));
      V = Vecteur(configFile.get<double>("vx0"),configFile.get<double>("vy0"),configFile.get<double>("vz0"));
      m = configFile.get<double>("m");
      R = configFile.get<double>("R");
      Cx= configFile.get<double>("Cx");
      corps=Corps(nom,P,V,m,R,Cx);

      systeme.add(corps);
    }
  }

  ~Exercice4()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.;
    last = 0;
    int i(0);
    cerr << "####################" << endl;
    printOut(true);
    while( t < tFin-0.5*dt )
    {
      systeme.evolue(t,dt,adaptatif,precision);
      printOut(false);
      if(i+1<t/tFin*20){
        i++;
        cerr << "#";
      }
    }
    cerr << "#" << endl;
    printOut(true);
  };

};


int main(int argc, char* argv[])
{ try{
  Exercice4 engine(argc, argv);
  engine.run();
}catch(Erreur err){
  err.affiche("Exercice4.cc");
}
  return 0;
}
