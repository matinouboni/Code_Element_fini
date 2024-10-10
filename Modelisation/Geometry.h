#ifndef _GEOMETRY_H

#include<vector>
#include<string>
#include "Dense"
#include "Sparse"
#include "mesh.h"

class Geometry
{
private:
  // Coordonnées des sommets du maillage
  Eigen::Matrix<double, Eigen::Dynamic, 2> _vertices;
  // Liste de tous les triangles et leur référence
  Eigen::Matrix<int, Eigen::Dynamic, 4> _triangles;
  // Liste de toutes les arêtes avec conditions de Neumann, leur référence,
  // leur triangle d'appartenance et le numéro de l'arête
  Eigen::Matrix<int, Eigen::Dynamic, 5> _edgesNeumann;

public:
  Geometry(Mesh* mesh);

  // Element de reference 
  // construction des fonction de base 
  // (phihat0, phihat1, phihat2)
  double phihat(int hati, Eigen::Vector2d);

  // Construction des 3 gradients des fonctions de base
  // (gradphihat0, gradphihat1, gradphihat2)
  // (indépendant du vecteur hataplpha)
  Eigen::Vector2d gradphihat(int hati);

  // Construction de la fonction de transformation
  // du triangle de référence hatalpha en alpha
  Eigen::Vector2d Falpha(int, Eigen::Vector2d);

  // Construction de la jacobienne de la fonction de transformation
  // du triangle de référence hatalpha en alpha
  Eigen::Matrix2d JFalpha(int);

  // Valeur absolue du determinant de la Jacobienne 
  double absoludetJ(int);

  // Construction de la fonction de transformation
  // de l'arête de référence hatE en E
  Eigen::Vector2d FE(int, Eigen::Vector2d);

  // Construction de la jacobienne de la fonction de transformation
  // de l'arête de référence hatE en E
  // (indépendant du vecteur hatX)
  Eigen::Matrix2d JFE(int);

  // Construction de la mesure associée à la transformation FE pour
  // le changement de variable de l'intégrale d Gamma_X = measK d Gamma_hatX
  double measE(int);




  // Construction de la matrice BoF_alpha (elementaire )
  Eigen::MatrixXd BoF_alpha(int);

  // Construction de la matrice D (en foncton de E et de nu)
  // de la loi de comportement 
  Eigen::MatrixXd coef_elastique(double, double);

  // calcul de la matrice produit 
  Eigen::MatrixXd Elementaire(int, double, double );

  // estimation des integrale 

  // Points et poids de quadrature de la formule du milieu (intégration 2D)
  void quadraturePointsAndWeightsMidpointFormula(Eigen::VectorXd& weights,
                                                 Eigen::Matrix<double, Eigen::Dynamic, 2>& points);

  // Points et poids de quadrature de la formule de Simpson (intégration 1D sur une arête)
  void quadraturePointsAndWeightsSimpsonFormula (Eigen::VectorXd& weights,
                                                 Eigen::Matrix<double, Eigen::Dynamic, 2>& points, int numedge);



  
};


#define _GEOMETRY_H
#endif
