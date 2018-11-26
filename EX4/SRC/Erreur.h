#pragma once
#include <string>

// classe Erreur qui aide la gestion des erreurs
class Erreur{
  public:
    //constructeur
    Erreur(std::string type="non spécifié", std::string fonction="non spécifiée", std::string dscrpt="Une erreur est survenue");
    // permet l'affichage de l'erreur
    void affiche(std::string program) const;
    //manipulateurs
    void add_fct(std::string newFct);
    void set_fct(std::string newFct);
    void set_dscrpt(std::string dscrpt);
  private:
    // attributs
    std::string type; // type d'erreur
	  std::string fct; // fonction dans laquelle l'erreur s'est produite
	  std::string description; // description de l'erreur
};
