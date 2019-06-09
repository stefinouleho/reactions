#ifndef NOM_D2
#define NOM_D2

#include "structure.h"
#include "graphe_cycles.h"
#include "helpers/graph.h"
#include "helpers/clique.h"

int valeur_absolue(int a);
GRAPHE_CYCLE ajouter_un_sommet(int a , int b, int taille ,GRAPHE_CYCLE c);
int type_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b);
int poids_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b);
GRAPHE_CYCLE ajouter_une_arete_graphe(int a , int b, int arete ,GRAPHE_CYCLE c);
GRAPHE_CYCLE construction_graphe_produit(GRAPHE_CYCLE a, GRAPHE_CYCLE b);
int ** construction_matrice_produit(GRAPHE_CYCLE p);
int somme_sommets_aretes( int *clique, GRAPHE_CYCLE produit, int **depart, int sommets);
void calcul_clique(int **matrice,int sommets,int *dans_clique,int taille_clique,int *candidat,int taille_candidat,GRAPHE_CYCLE produit, int **depart);
void la_clique_max( int **matrice, int sommets,int **depart,GRAPHE_CYCLE produit);
float similarite(GRAPHE_CYCLE a,GRAPHE_CYCLE b);
GRAPHE_CYCLE construction_commun_max(GRAPHE_CYCLE a,GRAPHE_CYCLE p);
void liberation_matrice(int **matrice, int taille);
graph graphe_produit_cycles(GRAPHE_CYCLE a,GRAPHE_CYCLE b);
struct couple *construction_couples_cycles(GRAPHE_CYCLE a,GRAPHE_CYCLE b,int taille);
struct type_arete ** construction_matrice_graphe_cycles(GRAPHE_CYCLE a);
void liberer_type_arete( struct type_arete **m, GRAPHE_CYCLE a);
int*  graphe_g12_cycles(graph g12, int* clique_max,GRAPHE_CYCLE a,GRAPHE_CYCLE b);

#endif