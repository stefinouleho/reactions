#ifndef NOM_D3
#define NOM_D3

#include "structure.h"

void init_atom_num ();
struct molecule * lecture_fichier_chebi();
struct molecule lire_molecule_sdf(FILE *F);
double chrono();
int lire_num_atome(FILE *F);
int valeur_char (FILE *F) ;
void ligne_suivante(FILE *F);
int lire_entier_3 (FILE * F);
struct liaison lire_liaison(FILE *F);
int lire_chebi_id(FILE *F);
void lire_chebi_name(FILE *F, struct molecule *M);
void lire_fin_molecule(FILE *F);
void lire_fin_molecule_2(FILE *F);
void trouver_la_fin_de_M(FILE *F);
void liberer_molecule(struct molecule m);
int atom_num (char *name);
int * lecture_liste_molecules();
int position_M( int g1_chebi,struct molecule *M);
void sauvegarde_compteur(int i , int j);

#endif
