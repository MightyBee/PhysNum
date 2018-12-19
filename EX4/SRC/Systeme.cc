#include <iostream>
#include <initializer_list>
#include <vector>
#include <cmath>
#include <string>
#include "Systeme.h"
using namespace std;


Corps::Corps(string nom, Vecteur P, Vecteur V, double m, double R, double Cx)
     : nom(nom), m(m), R(R), Cx(Cx) {y.P=P; y.V=V;}

ostream& Corps::affiche(ostream& sortie) const{
  sortie << y.P << " " << y.V;
  return sortie;
}

string Corps::getNom() const{
  return nom;
}

double Corps::getM() const{
  return m;
}

double Corps::getR() const{
  return R;
}

double Corps::getCx() const{
  return Cx;
}

VectY Corps::getY() const{
  return y;
}

void Corps::setY(const VectY& yRestaure){
  y=yRestaure;
}

ostream& operator<<(ostream& sortie, const Corps& corps){
  return corps.affiche(sortie);
}


Systeme::Systeme(double G, double rho0, double lambda)
       : G(G), rho0(rho0), lambda(lambda), liste(0) {}

ostream& Systeme::affiche(std::ostream& sortie) const{
  for(const auto& corps : liste){
    sortie << corps << " ";
  }
  return sortie;
}

void Systeme::add(const Corps& corps){
  liste.push_back(corps);
}

CollectionY Systeme::getY() const{
  CollectionY retour;
  for(const auto& corps : liste){
    retour.push_back(corps.getY());
  }
  return retour;
}

void Systeme::setY(const CollectionY& collY){
  for(size_t i(0); i<liste.size();i++){
    liste[i].setY(collY[i]);
  }
}

double Systeme::rho(double r) const{
  return rho0*exp((-(r-liste[0].getR()))/lambda);
}

CollectionY Systeme::f(const CollectionY& collY, double t) const{
  CollectionY retour(collY.taille());
  Vecteur FG;
  for(size_t i(0);i<collY.taille();i++){
    retour[i].P=collY[i].V;
    for(size_t k(i+1);k<collY.taille();k++){
      if(i!=k){
        FG=-G*(liste[i].getM()*liste[k].getM())/(collY[i].P-collY[k].P).norme2()*~(collY[i].P-collY[k].P);
        retour[i].V+=FG/liste[i].getM();
        retour[k].V-=FG/liste[k].getM();
      }
    }
    if(i!=0 and rho0!=0.0 and liste[i].getCx()!=0.0){
      retour[i].V+=-(1/liste[i].getM())*0.5*rho((collY[i].P-collY[0].P).norme())*M_PI*pow(liste[i].getR(),2)*liste[i].getCx()*(collY[i].V-collY[0].V).norme()*(collY[i].V-collY[0].V);
    }
  }
  return retour;
}


CollectionY Systeme::evolue(const CollectionY& y, double t, double dt){
  CollectionY k1=dt*f(y       ,t       );
  CollectionY k2=dt*f(y+0.5*k1,t+0.5*dt);
  CollectionY k3=dt*f(y+0.5*k2,t+0.5*dt);
  CollectionY k4=dt*f(y+    k3,t+    dt);
  CollectionY yNew=y+(k1+2*k2+2*k3+k4)/6;
  return yNew;
}

void Systeme::evolue(double& t, double& dt, bool adaptatif, double precision){
  if(adaptatif){
    double d(0.0);
    int i(0);
    CollectionY yOld=getY();
    CollectionY y2;
    do{
      if(i!=0){dt*=0.99*pow(precision/d,0.2);}
      y2=evolue(evolue(yOld,t,0.5*dt),t+0.5*dt,0.5*dt);
      d=(y2-evolue(yOld,t,dt)).norme();
      i++;
    }while(d>precision);
    setY(y2);
    t+=dt;
    dt*=pow(precision/d,0.2);
  }else{
    setY(evolue(getY(),t,dt));
    t+=dt;
  }
}

ostream& operator<<(ostream& sortie, const Systeme& syst){
  return syst.affiche(sortie);
}
