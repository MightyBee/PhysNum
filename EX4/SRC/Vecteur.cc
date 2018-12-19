#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <initializer_list>
#include "Vecteur.h"
using namespace std;


/*##############################################################################
###                                                                          ###
###                      METHODES DE LA CLASSE Vecteur                       ###
###                                                                          ###
##############################################################################*/

//#############################  constructeurs  ##############################//
// construit un vecteur nul de dimension n et fait office de constructeur par defaut
Vecteur::Vecteur(const double& val) : coord(3,val) {}

// construit un vecteur à partir d'une liste de double
Vecteur::Vecteur(const double& x, const double& y, const double& z) : coord({x,y,z}) {}


//##########################  opérateurs internes  ###########################//
// opérateurs de comparaison //
bool Vecteur::operator==(const Vecteur& v2) const{
	if(coord.size()!=v2.coord.size()){
		return false;
	} else {
		for(size_t i(0);i<coord.size();i++){
			if(coord[i]!=v2.coord[i]){return false;}
		}
		return true;
	}
}

bool Vecteur::operator!=(const Vecteur& v2) const{
	return not operator==(v2);
}


// retourne produit scalaire du vecteur courant avec un autre vecteur //
double Vecteur::operator*(const Vecteur& v2) const{
	if(coord.size()==v2.coord.size()){
		double retour(0.0);
		for(size_t i(0);i<coord.size();i++){
			retour+=coord[i]*v2.coord[i];
		}
		return retour;
	} else {
		Erreur err("dimension", "Vecteur::operator*(const Vecteur&)",
							 "Produit scalaire de deux vecteurs de dimension différente ("+to_string(coord.size())+" et "+to_string(v2.coord.size())+").");
		throw err;
	}
}

// addition d'un vecteur au vecteur courant //
Vecteur& Vecteur::operator+=(const Vecteur& v2){
	if(coord.size()==v2.coord.size()){
		for(size_t i(0);i<coord.size();i++){
			coord[i]+=v2.coord[i];
		}
		return *this;
	} else {
		Erreur err("dimension", "Vecteur::operator+=(const Vecteur&)",
							 "Addition/soustraction de deux vecteurs de dimension différente ("+to_string(coord.size())+" et "+to_string(v2.coord.size())+").");
		throw err;
	}
}

// soustraction d'un vecteur au vecteur courant //
Vecteur& Vecteur::operator-=(const Vecteur& v2){
	try{
		operator+=(-v2);
		return *this;
	}catch(Erreur err){
		err.set_fct("Vecteur::operator-=(const Vecteur&)");
		throw err;
	}
}

// multiplication du vecteur courant par un scalaire //
Vecteur& Vecteur::operator*=(const double& lambda){
	for(auto& el : coord){el*=lambda;}
	return *this;
}

// division du vecteur courant par un scalaire //
Vecteur& Vecteur::operator/=(const double& lambda){
	if(lambda!=0){
		operator*=(1.0/lambda);
		return *this;
	} else {
		Erreur err("division par 0", "Vecteur::operator/=(const double&)", "Division d'un vecteur par zéro.");
		throw err;
	}
}


//##############################  accesseurs  ################################//
const double& Vecteur::operator[](unsigned int i) const{
	if(i<coord.size()){
		return coord[i];
	}else{
		Erreur err("dimension", "Vecteur::operator[](unsigned int) const",
							 "L'indice de position fourni en argument ("+to_string(i)+") n'est pas valable (attendu : entier entre 0 et "+to_string(coord.size()-1)+", dim(Vecteur)="+to_string(coord.size())+").");
		throw err;
	}
}

// retourne la dimension du vecteur //
size_t Vecteur::taille() const{
	return coord.size();
}


//#############################  manipulateurs  ##############################//
// retourne une référence sur la (i+1)-ième coordonnee du vecteur, permet de la modifier //
double& Vecteur::operator[](unsigned int i){
	if(i>=0 and i<coord.size()){
		return coord[i];
	}else{
		Erreur err("dimension", "Vecteur::operator[](unsigned int)",
							 "L'indice de position fourni en argument ("+to_string(i)+") n'est pas valable (attendu : entier entre 0 et "+to_string(coord.size()-1)+", dim(Vecteur)="+to_string(coord.size())+").");
		throw err;
	}
}


//###########################  autres opérations  ############################//
// affichage des coordonnees du vecteur //
ostream& Vecteur::affiche(ostream& sortie)const{
		sortie << x() << " " << y() << " " << z() ;
	return sortie;
}

// retourne la norme du vecteur courant //
double Vecteur::norme() const{
	return sqrt(norme2());
}

// retourne la norme au carré du vecteur courant : c'est <v,v> //
double Vecteur::norme2() const{
	return (*this)*(*this); // on utilise la surcharge de l'operateur * (produit sclaire)
}


// retourne la projection du vecteur sur le plan XY
Vecteur Vecteur::projXY() const{
	return Vecteur(x(),y(),0);
}


// retourne l'angle en valeur absolue entre deux vecteur //
double Vecteur::angle(Vecteur const& v2) const{
	double angle(0);
	Vecteur V1(~(*this)); // on les normalise pour que <·,·>=cos(angle)
	Vecteur V2(~v2);
	if(V1==-V2){angle=-M_PI;}
  else if(V1!=V2){
    angle=acos(V1*V2);
  }
	return angle;
}


/*##############################################################################
###                                                                          ###
###                           OPÉRATEURS EXTERNES                            ###
###                                                                          ###
##############################################################################*/

// permet l'affichage standard : sortie << vecteur; //
ostream& operator<<(ostream& sortie, const Vecteur& v){
	return v.affiche(sortie);
}

// somme de deux vecteurs //
const Vecteur operator+(Vecteur v1, const Vecteur& v2){
	try{
		v1+=v2;
		return v1;
	}catch(Erreur err){
		err.set_fct("operator+(Vecteur, const Vecteur&)");
		throw err;
	}
}

// différence de deux vecteurs //
const Vecteur operator-(Vecteur v1, const Vecteur& v2){
	try{
		v1-=v2;
		return v1;
	}catch(Erreur err){
		err.set_fct("operator-(Vecteur, const Vecteur&)");
		throw err;
	}
}

// opposé d'un vecteur //
const Vecteur operator-(Vecteur v){
	v*=-1;
	return v;
}

// multiplication d'un vecteur par un scalaire, cas scal*vect //
const Vecteur operator*(const double& lambda, Vecteur v){
	v*=lambda;
	return v;
}

// multiplication d'un vecteur par un scalaire, cas vect*scal //
const Vecteur operator*(const Vecteur& v, const double& lambda){
	return lambda*v;
}


// division d'un vecteur par un scalaire //
const Vecteur operator/(Vecteur v, double lambda){
	try{
		v/=lambda;
		return v;
	}catch(Erreur err){
		err.set_fct("operator/(Vecteur, double)");
		throw err;
	}
}


// retourne le vecteur unitaire : v/||v||
const Vecteur operator~(Vecteur v){
	try{
		return v/v.norme();
	}catch(Erreur err){ // on catch l'erreur division par zéro pour la relancer
		err.set_fct("operator~(Vecteur), i.e. vecteur unitaire");
		err.set_dscrpt("Le vecteur nul n'a pas de vecteur unitaire correspondant.");
		throw err;
	}
}


// retourne le produit vectoriel du vecteur courant avec un autre vecteur 3D //
const Vecteur operator^(const Vecteur& v1, const Vecteur& v2){
	return Vecteur(v1.y()*v2.z()-v1.z()*v2.y(),
								 v1.z()*v2.x()-v1.x()*v2.z(),
								 v1.x()*v2.y()-v1.y()*v2.x());
}




CollectionY::CollectionY(const double& n) : collection(n) {
	for(auto& y : collection){
		y.P=Vecteur(0);
		y.V=Vecteur(0);
	}
}

// opérateurs de comparaison //
bool CollectionY::operator==(const CollectionY& collY2) const{
	if(collection.size()!=collY2.collection.size()){
		return false;
	} else {
		for(size_t i(0);i<collection.size();i++){
			if(collection[i].P!=collY2.collection[i].P or collection[i].V!=collY2.collection[i].V){
        return false;
      }
		}
		return true;
	}
}

bool CollectionY::operator!=(const CollectionY& collY2) const{
	return not operator==(collY2);
}


// addition d'un vecteur au vecteur courant //
CollectionY& CollectionY::operator+=(const CollectionY& collY2){
	if(collection.size()==collY2.collection.size()){
		for(size_t i(0);i<collection.size();i++){
			collection[i].P+=collY2.collection[i].P;
      collection[i].V+=collY2.collection[i].V;
		}
		return *this;
	} else {
		Erreur err("dimension", "CollectionY::operator+=(const CollectionY&)",
							 "Addition/soustraction de deux vecteurs de dimension différente ("+to_string(collection.size())+" et "+to_string(collY2.collection.size())+").");
		throw err;
	}
}

// soustraction d'un vecteur au vecteur courant //
CollectionY& CollectionY::operator-=(const CollectionY& collY2){
	try{
		operator+=(-collY2);
		return *this;
	}catch(Erreur err){
		err.set_fct("CollectionY::operator-=(const CollectionY&)");
		throw err;
	}
}

// multiplication du vecteur courant par un scalaire //
CollectionY& CollectionY::operator*=(const double& lambda){
	for(auto& el : collection){
    el.P*=lambda;
    el.V*=lambda;
  }
	return *this;
}

// division du vecteur courant par un scalaire //
CollectionY& CollectionY::operator/=(const double& lambda){
	if(lambda!=0){
		operator*=(1.0/lambda);
		return *this;
	} else {
		Erreur err("division par 0", "CollectionY::operator/=(const double&)", "Division d'un vecteur par zéro.");
		throw err;
	}
}


//##############################  accesseurs  ################################//
const VectY& CollectionY::operator[](unsigned int i) const{
	if(i<collection.size()){
		return collection[i];
	}else{
		Erreur err("dimension", "CollectionY::operator[](unsigned int) const",
							 "L'indice de position fourni en argument ("+to_string(i)+") n'est pas valable (attendu : entier entre 0 et "+to_string(collection.size()-1)+", dim(CollectionY)="+to_string(collection.size())+").");
		throw err;
	}
}

// retourne la dimension du vecteur //
size_t CollectionY::taille() const{
	return collection.size();
}

void CollectionY::push_back(const VectY& y){
  collection.push_back(y);
}

//#############################  manipulateurs  ##############################//
// retourne une référence sur la (i+1)-ième coordonnee du vecteur, permet de la modifier //
VectY& CollectionY::operator[](unsigned int i){
	if(i>=0 and i<collection.size()){
		return collection[i];
	}else{
		Erreur err("dimension", "CollectionY::operator[](unsigned int)",
							 "L'indice de position fourni en argument ("+to_string(i)+") n'est pas valable (attendu : entier entre 0 et "+to_string(collection.size()-1)+", dim(CollectionY)="+to_string(collection.size())+").");
		throw err;
	}
}


//###########################  autres opérations  ############################//
// retourne la norme du vecteur courant //
double CollectionY::norme() const{
	return sqrt(norme2());
}

// retourne la norme au carré du vecteur courant : c'est <v,v> //
double CollectionY::norme2() const{
  double retour;
  for(const auto& y : collection){
    retour+=y.P.norme2()+y.V.norme2();
  }
	return retour;
}

// affichage des coordonnees du vecteur //
ostream& CollectionY::affiche(ostream& sortie)const{
	//sortie << "( ";
	for(auto el : collection){
		sortie << el.P << " " << el.V ;
	}
	//sortie << ")";
	return sortie;
}



/*##############################################################################
###                                                                          ###
###                           OPÉRATEURS EXTERNES                            ###
###                                                                          ###
##############################################################################*/

// permet l'affichage standard : sortie << vecteur; //
ostream& operator<<(ostream& sortie, const CollectionY& collY){
	return collY.affiche(sortie);
}

// somme de deux vecteurs //
const CollectionY operator+(CollectionY collY1, const CollectionY& collY2){
	try{
		collY1+=collY2;
		return collY1;
	}catch(Erreur err){
		err.set_fct("operator+(CollectionY, const CollectionY&)");
		throw err;
	}
}

// différence de deux vecteurs //
const CollectionY operator-(CollectionY collY1, const CollectionY& collY2){
	try{
		collY1-=collY2;
		return collY1;
	}catch(Erreur err){
		err.set_fct("operator-(CollectionY, const CollectionY&)");
		throw err;
	}
}

// opposé d'un vecteur //
const CollectionY operator-(CollectionY collY){
	collY*=-1;
	return collY;
}

// multiplication d'un vecteur par un scalaire, cas scal*vect //
const CollectionY operator*(const double& lambda, CollectionY collY){
	collY*=lambda;
	return collY;
}

// multiplication d'un vecteur par un scalaire, cas vect*scal //
const CollectionY operator*(const CollectionY& collY, const double& lambda){
	return lambda*collY;
}


// division d'un vecteur par un scalaire //
const CollectionY operator/(CollectionY collY, double lambda){
	try{
		collY/=lambda;
		return collY;
	}catch(Erreur err){
		err.set_fct("operator/(CollectionY, double)");
		throw err;
	}
}
