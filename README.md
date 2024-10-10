# Modélisation à partir d’images 2D du comportement mécanique de composites : étude à l'échelle des torons

## Description du projet

Ce projet a pour objectif de modéliser le comportement mécanique des composites à l'échelle des torons, en utilisant une approche par éléments finis. La méthode repose sur l'utilisation d'un élément de référence pour simplifier et systématiser les calculs, notamment pour la construction des matrices de rigidité et des termes sources.

### Passage par élément de référence

Pour simplifier les calculs, nous utilisons un **élément de référence** : un triangle avec des sommets définis en \((0,0)\), \((1,0)\) et \((0,1)\). Tout triangle du maillage est alors une image affine de cet élément de référence. Cela permet de simplifier les intégrations numériques et les transformations géométriques, facilitant ainsi l'assemblage du système global.

### Formulation discrète

La formulation discrète des équations mécaniques est basée sur l'approximation par fonctions de base définies sur les éléments de référence. Les efforts internes sont calculés à partir des déformations et des contraintes, à travers des relations de comportement du matériau, et sont discrétisés en un système d'équations résolu numériquement.

Le projet intègre également des approximations numériques, telles que l'utilisation de quadratures de Gauss pour résoudre les intégrales dans les éléments finis.

### Principales étapes du code :

1. **Lecture des données d'entrée** : Le fichier `data.txt` contient les paramètres matériels et les conditions aux limites.
2. **Maillage** : Le maillage est généré à partir des fichiers contenus dans le dossier `Mesh/`.
3. **Calcul géométrique** : Les matrices de transformation et de déformation sont calculées sur chaque élément du maillage.
4. **Assemblage du système global** : En tenant compte des conditions aux limites, le système d'équations global est construit.
5. **Résolution** : Le solveur utilise les méthodes définies dans `Solvuer.cpp` pour résoudre le système linéaire.
6. **Sortie des résultats** : Les résultats (déplacements, déformations, contraintes) sont exportés dans `Resultats.vtk` pour visualisation.

## Structure du projet

- `main.cc` : Le point d'entrée principal du programme.
- `Mesh/` : Contient les fichiers de maillage utilisés pour les simulations.
- `Geometry.cpp`, `Geometry.h` : Contiennent les calculs géométriques pour les éléments finis.
- `Assemblage.cpp`, `Assemblage.h` : Gèrent l'assemblage des matrices de rigidité et du second membre.
- `Solvuer.cpp`, `Solvuer.h` : Contiennent le solveur du système linéaire.
- `mesh.cpp`, `mesh.h` : Gèrent les données de maillage.
- `data.txt` : Données d'entrée pour les simulations.


## Remarque : la bibliothèque Eigen est nécessaire pour la compilation

## Compilation et Exécution

### Compilation

Un fichier **Makefile** est fourni pour compiler l'ensemble du projet. Pour compiler le projet, exécutez la commande suivante dans un terminal à la racine du projet :

```bash
make
 ./run data.txt
 
 Vous pouvez modifier le fichier data.txt en fonction de vos paramètres d'entrée spécifiques.
 

