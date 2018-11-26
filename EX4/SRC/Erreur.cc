#include "Erreur.h"
#include <iostream>
#include <string>
using namespace std;

// constructeur
Erreur::Erreur(std::string type, std::string fonction, std::string dscrpt)
              : type(type), fct(fonction), description(dscrpt) {}

// affichage du message d'erreur
void Erreur::affiche(string program) const {
  cerr << endl << endl << "### ERREUR FATALE ###" << endl;
  cerr << "Type d'erreur : " << type << endl;
  cerr << "Dans fichier  : "+program << endl;
  cerr << "Dans fonction : " << fct << endl;
  cerr << "Description   : " << description << endl;
  cerr << "### FIN DU PROGRAMME ###" << endl << endl;
}

// rajout d'une fonction: losqu'une erruer est produite par une fonction appelée par une autre fonction
void Erreur::add_fct(string newFct){
  fct+=" appelée par "+newFct;
}

void Erreur::set_fct(std::string newFct){
  fct=newFct;
}

void Erreur::set_dscrpt(std::string dscrpt){
  description=dscrpt;
}
