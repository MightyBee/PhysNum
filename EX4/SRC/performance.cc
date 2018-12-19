#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <array>
#include <string>
#include <initializer_list>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template
#include <valarray>
using namespace std;


double norme(const valarray<double>& y){
  return sqrt(pow(y,2.0).sum());
}


class Exercice4
{

private:
  double t, dt, precision, tFin;
  double G, rho0, lambda;
  bool adaptatif;
  array<double,3> m, R, Cx;
  valarray<double> y;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(const bool& force)
  {
    if((!force && last>=sampling) || (force && last!=1)){
      *outputFile << t << " " << accAstronautes() << " " << emec() << " " << PT() << " ";
      *outputFile << y[0] << " " << y[1] << " 0 " << y[2]  << " " << y[3]  << " 0 ";
      *outputFile << y[4] << " " << y[5] << " 0 " << y[6]  << " " << y[7]  << " 0 ";
      *outputFile << y[8] << " " << y[9] << " 0 " << y[10] << " " << y[11] << " 0 ";
      *outputFile << endl;
      last = 1;
    }
    else{
      last++;
    }
  }

  double accAstronautes() const{
    valarray<double> acc(2);
    valarray<double> r0=y[slice(0,2,1)];
    valarray<double> r1=y[slice(4,2,1)];
    valarray<double> r2=y[slice(8,2,1)];
    acc=-G*(m[0]/pow(pow(r2-r0,2.0).sum(),1.5)*(r2-r0)+m[1]/pow(pow(r2-r1,2.0).sum(),1.5)*(r2-r1));
    if(rho0!=0.0){
      acc += -(0.5/m[2])*rho(norme(r2-r0))*M_PI*R[2]*R[2]*Cx[2]*norme(y[slice(10,2,1)]-y[slice(2,2,1)])*(y[slice(10,2,1)]-y[slice(2,2,1)]);
    }
    return norme(acc);
  }

  double emec() const{
    double retour(0.0);
    valarray<double> r0=y[slice(0,2,1)];
    valarray<double> r1=y[slice(4,2,1)];
    valarray<double> r2=y[slice(8,2,1)];
    retour-=G*m[0]*m[1]/norme(r0-r1);
    retour-=G*m[0]*m[2]/norme(r0-r2);
    retour-=G*m[1]*m[2]/norme(r1-r2);
    return retour;
  }

  double rho(const double& r) const{
    return rho0*exp(-(r-R[0])/lambda);
  }

double PT() const{
  valarray<double> r0=y[slice(0 ,2,1)];
  valarray<double> v0=y[slice(2 ,2,1)];
  valarray<double> r2=y[slice(8 ,2,1)];
  valarray<double> v2=y[slice(10,2,1)];
  if(rho0!=0.0){
    return (v2*(-0.5*rho(norme(r2-r0))*M_PI*R[2]*R[2]*Cx[2]*norme(v2-v0)*(v2-v0))).sum();
}else{
  return 0;
}
}

  valarray<double> f(const valarray<double>& y, const double& t) const{
    valarray<double> retour(12);
    retour[slice(0,3,4)]=y[slice(2,3,4)];
    retour[slice(1,3,4)]=y[slice(3,3,4)];
    valarray<double> r0=y[slice(0,2,1)];
    valarray<double> r1=y[slice(4,2,1)];
    valarray<double> r2=y[slice(8,2,1)];
    retour[slice(2,2,1)]  = -G*(m[1]/pow(pow(r0-r1,2).sum(),1.5)*(r0-r1)+m[2]/pow(pow(r0-r2,2).sum(),1.5)*(r0-r2));
    retour[slice(6,2,1)]  = -G*(m[0]/pow(pow(r1-r0,2).sum(),1.5)*(r1-r0)+m[2]/pow(pow(r1-r2,2).sum(),1.5)*(r1-r2));
    retour[slice(10,2,1)] = -G*(m[0]/pow(pow(r2-r0,2).sum(),1.5)*(r2-r0)+m[1]/pow(pow(r2-r1,2).sum(),1.5)*(r2-r1));
    if(rho0!=0.0){
      retour[slice(10,2,1)]+= -(0.5/m[2])*rho(norme(r2-r0))*M_PI*R[2]*R[2]*Cx[2]*norme(y[slice(10,2,1)]-y[slice(2,2,1)])*(y[slice(10,2,1)]-y[slice(2,2,1)]);
    }
    return retour;
  }

  valarray<double> one_step(const valarray<double>& y, const double& t, const double& dt) const {
    valarray<double> k1=dt*f(y       ,t       );
    valarray<double> k2=dt*f(y+0.5*k1,t+0.5*dt);
    valarray<double> k3=dt*f(y+0.5*k2,t+0.5*dt);
    valarray<double> k4=dt*f(y+    k3,t+    dt);
    return y+(k1+2.0*k2+2.0*k3+k4)/6;
  }

  void evolue(){
    if(t+dt>tFin){
      dt=tFin-t;
    }
    if(adaptatif){
      double d(0.0), fact(1.0);
      int i(0);
      valarray<double> y2(12);
      do{
        if(i!=0){dt*=0.97*pow(precision*fact/d,0.2);}
        y2=one_step(one_step(y,t,0.5*dt),t+0.5*dt,0.5*dt);
        //fact/=norme(fact*y2);
        d=norme(fact*y2-fact*one_step(y,t,dt));
        i++;
      }while(d>precision*fact);
      y=y2;
      t+=dt;
      dt*=pow(precision*fact/d,0.2);
    }else{
      y=one_step(y,t,dt);
      t+=dt;
    }
  }



public:

  Exercice4(int argc, char* argv[]) : y(12)
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
    G         = configFile.get<double>("G");
    rho0      = configFile.get<double>("rho0");
    lambda    = configFile.get<double>("lambda");
    sampling  = configFile.get<int>("sampling");
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    // Ajout des corps
    for(size_t i(0); i<3; i++){
      inputPath0="configuration"+to_string(i+1)+".in"; // Fichier d'input par defaut
      ConfigFile configFile("config/"+inputPath0); // Les parametres sont lus et stockes dans une "map" de strings.
      for(size_t k(0);k<inputPath.size();k++){
        if(inputPath[k]==inputPath0){
          for(const auto& str : param[k]){
            configFile.process(str);
          }
        }
      }
      y[4*i]   = configFile.get<double>("x0");
      y[4*i+1] = configFile.get<double>("y0");
      y[4*i+2] = configFile.get<double>("vx0");
      y[4*i+3] = configFile.get<double>("vy0");
      m[i]     = configFile.get<double>("m");
      R[i]     = configFile.get<double>("R");
      Cx[i]    =configFile.get<double>("Cx");
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
    //int i(0);
    //cerr << "####################" << endl;
    printOut(true);
    while( t < tFin and norme(valarray<double>(y[slice(8,2,1)])-valarray<double>(y[slice(0,2,1)])) > R[0] and norme(valarray<double>(y[slice(8,2,1)])-valarray<double>(y[slice(4,2,1)])) >= R[1])
    {
      evolue();
      printOut(false);
      //if(i+1<t/tFin*20){
      //  i++;
      //  cerr << "#";
      //}
    }
    //cerr << "#" << endl;
    printOut(true);
  };

};


int main(int argc, char* argv[])
{
  Exercice4 engine(argc, argv);
  engine.run();
  return 0;
}
