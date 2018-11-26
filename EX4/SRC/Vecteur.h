#pragma once
#include "Erreur.h"
#include <vector>
#include <string>
#include <cmath>
#include <initializer_list>


/*##############################################################################
###                                                                          ###
###                             CLASSE Vecteur                               ###
###                                                                          ###
##############################################################################*/

class Vecteur {
	public:
		//constructeurs
		explicit Vecteur(const double& val=0); // construit un vecteur de dimension n et fait office de constructeur par defaut
		Vecteur(const double& x, const double& y, const double& z); // construit un vecteur à partir d'une liste de double

		//operateurs internes
		bool operator==(const Vecteur& v2) const; // retourne true si les vecteurs sont les memes, false sinon
		bool operator!=(const Vecteur& v2) const; // retourne true si les vecteurs sont differents, flase sinon

		double operator*(const Vecteur& v2) const; // retourne le produit scalaire de l'instance courante avec un autre vecteur

		Vecteur& operator+=(const Vecteur& v2); // additionne un vecteur à l'instance courante
		Vecteur& operator-=(const Vecteur& v2); // soustrait un vecteur à l'instance courante
		Vecteur& operator*=(const double& lambda); // multiplie l'instance courante par un scalaire (double)
		Vecteur& operator/=(const double& lambda); // divise l'instance courante par un scalaire (double)

		//accesseurs
		const double& operator[](unsigned int i) const; // accesseur, ième coordonnée
		double x() const{return coord[0];}
		double y() const{return coord[1];}
		double z() const{return coord[2];}
		size_t taille() const; // retourne la dim du vecteur

		//manipulateurs
		double& operator[](unsigned int i); // accesseur et manipulateur, (i+1)ème coordonnée

		//autres operations
		double norme() const; // retourne la norme de l'instance courante
		double norme2() const; // retourne la norme au carré de l'instance courante
		Vecteur projXY() const; // retourne la projection du vecteur sur le plan XY
		double angle(Vecteur const& v2) const; // retourne l'angle entre un vecteur et l'instance courante
		std::ostream& affiche(std::ostream& sortie) const; // permet d'afficher le vecteur par composants sur un flot de sortie

	private:
		//attributs
		std::vector<double> coord; // tableau dynamique contenant les coordonees du vecteur
};


// operateurs externes
std::ostream& operator<<(std::ostream& sortie, const Vecteur& v); // permet l'affichage standard : sortie << vecteur;
const Vecteur operator+(Vecteur v1, const Vecteur& v2); // somme de deux vecteurs
const Vecteur operator-(Vecteur v1, const Vecteur& v2); // différence de deux vecteurs
const Vecteur operator-(Vecteur v); // opposé d'un vecteur
const Vecteur operator*(const double& lambda, Vecteur v); // multiplication d'un vecteur par un scalaire, cas scal*vect
const Vecteur operator*(const Vecteur& v, const double& lambda); // multiplication d'un vecteur par un scalaire, cas vect*scal
const Vecteur operator/(Vecteur v, double lambda); // division d'un vecteur par un scalaire
const Vecteur operator^(const Vecteur& v1, const Vecteur& v2); // retourne le produit vectoriel de l'instance courante avec un autre vecteur 3D
const Vecteur operator~(Vecteur v); // retourne le vecteur unitaire : v/||v||

struct VectY{
    Vecteur P;
    Vecteur V;
};

class CollectionY{
  public:
    //operateurs internes
		bool operator==(const CollectionY& collY2) const; // retourne true si les vecteurs sont les memes, false sinon
		bool operator!=(const CollectionY& collY2) const; // retourne true si les vecteurs sont differents, false sinon
		CollectionY& operator+=(const CollectionY& collY2); // additionne un vecteur à l'instance courante
		CollectionY& operator-=(const CollectionY& collY2); // soustrait un vecteur à l'instance courante
		CollectionY& operator*=(const double& lambda); // multiplie l'instance courante par un scalaire (double)
		CollectionY& operator/=(const double& lambda); // divise l'instance courante par un scalaire (double)
		//accesseurs
    size_t taille() const;
		const VectY& operator[](unsigned int i) const; // accesseur, ième coordonnée
    //manipulateurs
    void push_back(const VectY& y);
		VectY& operator[](unsigned int i); // accesseur et manipulateur, (i+1)ème coordonnée
    //autres operations
    double norme() const; // retourne la norme de l'instance courante
		double norme2() const; // retourne la norme au carré de l'instance courante
		std::ostream& affiche(std::ostream& sortie) const; // permet d'afficher le vecteur par composants sur un flot de sortie
  private:
    std::vector<VectY> collection;
};

// operateurs externes
std::ostream& operator<<(std::ostream& sortie, const CollectionY& v); // permet l'affichage standard : sortie << vecteur;
const CollectionY operator+(CollectionY collY1, const CollectionY& collY2); // somme de deux vecteurs
const CollectionY operator-(CollectionY collY1, const CollectionY& collY2); // différence de deux vecteurs
const CollectionY operator-(CollectionY collY); // opposé d'un vecteur
const CollectionY operator*(const double& lambda, CollectionY collY); // multiplication d'un vecteur par un scalaire, cas scal*vect
const CollectionY operator*(const CollectionY& collY, const double& lambda); // multiplication d'un vecteur par un scalaire, cas vect*scal
const CollectionY operator/(CollectionY collY, double lambda); // division d'un vecteur par un scalaire
