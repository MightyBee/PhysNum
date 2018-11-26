#pragma once
#include <initializer_list>
#include <vector>
#include <string>
#include "Vecteur.h"
#include "Erreur.h"

class Corps{
  public:
    Corps(){}
    Corps(std::string nom, Vecteur P, Vecteur V, double m, double R, double Cx);
    std::ostream& affiche(std::ostream& sortie) const;
    std::string getNom() const;
    double getM() const;
    double getR() const;
    double getCx() const;
    VectY getY() const;
    void setY(const VectY& yRestaure);
  private:
    std::string nom; // nom du corps
    double m;  // masse du corps
    double R;  // rayon du corps
    double Cx; // coef de trainee
    VectY y; // position et vitesse du corps
};
std::ostream& operator<<(std::ostream& sortie, const Corps& corps);

class Systeme{
  public:
    Systeme(){}
    Systeme(double G, double rho0, double lambda);
    std::ostream& affiche(std::ostream& sortie) const;
    void add(const Corps& corps);
    CollectionY getY() const;
    void setY(const CollectionY& s);
    double rho(double r) const;
    CollectionY f(const CollectionY& collY, double t) const;
    CollectionY evolue(double t, double dt);
    void evolue(double& t, double& dt, bool adaptatif, double precision);
  private:
    double G;            // constante de gravitation
    double rho0;         // densité de l'air au niveau de la mer
    double lambda;       // épaisseur caractéristique
    std::vector<Corps> liste; // liste des corps
};
std::ostream& operator<<(std::ostream& sortie, const Systeme& syst);
