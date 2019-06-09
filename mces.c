#include "mces.h"


int taille_clique_max = 0;
int *dans_clique_max= NULL;


int valeur_absolue(int a)
{
	if ( a > 0)
		return a;
	return -a;
}

GRAPHE_CYCLE ajouter_un_sommet(int a , int b, int taille ,GRAPHE_CYCLE c)
{
	if( c.liste_sommets == NULL)
	{
		c.liste_sommets = malloc((taille + 1)* sizeof(SOMMET));
	}
	else
	{
		c.liste_sommets = realloc(c.liste_sommets,(taille +1)* sizeof(SOMMET));

	}
	if(c.liste_sommets == NULL)
		probleme_memoire();
	c.liste_sommets[taille].id = a;
	c.liste_sommets[taille].poids = b;
	return c;
}

int type_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b)
{
	int res = -1;

	int i;
	for ( i = 0; i < c.nb_aretes;i++)
	{
		if((c.liste_aretes[i].id1 == a && c.liste_aretes[i].id2 == b)||(c.liste_aretes[i].id1 == b && c.liste_aretes[i].id2 == a))
		{
			res = c.liste_aretes[i].type;
			break;
		}	
	}
	return res;
}
int poids_arete_graphe_cycle(GRAPHE_CYCLE c,int a ,int b)
{
	int res = -1;

	int i;
	for ( i = 0; i < c.nb_aretes;i++)
	{
		if((c.liste_aretes[i].id1 == a && c.liste_aretes[i].id2 == b)||(c.liste_aretes[i].id1 == b && c.liste_aretes[i].id2 == a))
		{
			res = c.liste_aretes[i].poids;
			break;
		}	
	}
	return res;
}

GRAPHE_CYCLE ajouter_une_arete_graphe(int a , int b, int arete ,GRAPHE_CYCLE c)
{
	if( c.liste_aretes == NULL)
	{
		c.liste_aretes = malloc((arete + 1)* sizeof(ARETE));
	}
	else
	{
		c.liste_aretes= realloc(c.liste_aretes,(arete +1)* sizeof(ARETE));

	}
	if(c.liste_aretes == NULL)
		probleme_memoire();
	c.liste_aretes[arete].id1 = a;
	c.liste_aretes[arete].id2 = b;
	c.liste_aretes[arete].poids = 1;
	return c;
}

GRAPHE_CYCLE construction_graphe_produit(GRAPHE_CYCLE a, GRAPHE_CYCLE b)
{
	GRAPHE_CYCLE c;
	c.liste_aretes = NULL;
	c.liste_sommets = NULL;
	int taille = 0;
	//printf("nb cycles de a %d et nb cycles de b %d\n",a.nb_sommets,b.nb_sommets );
	int i,j;
	for ( i = 0; i < a.nb_sommets; i++)
	{
		//printf("%d :",a.liste_sommets[i].poids );
		for ( j = 0; j < b.nb_sommets; j++)
		{
			if( (float)(valeur_absolue(a.liste_sommets[i].poids - b.liste_sommets[j].poids)) <= 0.2 *min(a.liste_sommets[i].poids ,b.liste_sommets[j].poids))
			{
				c = ajouter_un_sommet(a.liste_sommets[i].id,b.liste_sommets[j].id,taille,c );
				taille++;
				//printf("%d ",b.liste_sommets[j].poids );
			}
		}
	}
	c.nb_sommets = taille;
	int aretes = 0;
	c.liste_aretes = NULL;
	int res1,res2,res3 = 0,res4 = 1;
	for ( i = 0; i < c.nb_sommets - 1; i++)
	{
		for ( j = i + 1; j < c.nb_sommets; j++)
		{
			res3 = 0;
			res4 = 1;
			//printf("le couplet(%d, %d) --- (%d, %d)\n",c.liste_sommets[i].id,c.liste_sommets[i].poids,c.liste_sommets[j].id,c.liste_sommets[j].poids );
			
			if(c.liste_sommets[i].id == c.liste_sommets[j].id)
			{
				res1 = -1;
			}
			else
			{
				res1 = type_arete_graphe_cycle(a,c.liste_sommets[i].id , c.liste_sommets[j].id);
				res3 = poids_arete_graphe_cycle(a,c.liste_sommets[i].id , c.liste_sommets[j].id);
			}	
			if(c.liste_sommets[i].poids == c.liste_sommets[j].poids)
			{
				res2 = -1;
			}
			else
			{
				res2 = type_arete_graphe_cycle(b,c.liste_sommets[i].poids , c.liste_sommets[j].poids);
				res4 = poids_arete_graphe_cycle(b,c.liste_sommets[i].poids , c.liste_sommets[j].poids);
			}	
			if( res1 * res2 >= 0)
			{
				if(res1 == res2 && res3 == res4)
				{
					c = ajouter_une_arete_graphe(i , j,aretes,c);
					aretes++;
					//printf("le couplet(%d, %d) --- (%d, %d) type = %d et poids %d\n",c.liste_sommets[i].id,c.liste_sommets[i].poids,c.liste_sommets[j].id,c.liste_sommets[j].poids,res1,res3 );
			
				}
				/*else 
				{
					if(res1 == res2 && (res1 == 2 || res1 == 3))
					{
						c = ajouter_une_arete_graphe(i , j,aretes,c);
						aretes++;
					}
				}*/
			}
		}
	}
	c.nb_aretes = aretes;
	//printf("taille du graphe produit %d et aretes %d \n", c.nb_sommets, c.nb_aretes);

	return c;
}

int ** construction_matrice_produit(GRAPHE_CYCLE p)
{
	int i,j;
	int **matrice = malloc(p.nb_sommets * sizeof(int *));
	if( matrice == NULL)
		probleme_memoire();
	for(i = 0; i < p.nb_sommets;i++)
	{
		matrice[i] = malloc(p.nb_sommets * sizeof(int));
	 	if(matrice[i] == NULL)
	 		probleme_memoire();
	}
	for(i = 0; i < p.nb_sommets;i++)
	{
		for(j = 0; j < p.nb_sommets;j++)
	 	{
	 		matrice[i][j] = 0;
	 	}
	}
	for(i = 0; i < p.nb_aretes;i++)
	{
	 	matrice[p.liste_aretes[i].id1][p.liste_aretes[i].id2] = 1;
	 	matrice[p.liste_aretes[i].id2][p.liste_aretes[i].id1] = 1;
	}
	 return matrice;
}


int somme_sommets_aretes( int *clique, GRAPHE_CYCLE produit, int **depart, int sommets)
{

	//printf("i am here \n");
	int somme = 0;
	int i,j;
	
	for( i = 0; i < sommets;i++)
	{
		somme += clique[i];
	}
	for( i = 0; i < sommets - 1 ; i++)
	{
		for( j = i+1; j < sommets; j++)
		{
			
			int id1,id2;
			id1 = produit.liste_sommets[i].id;
			id2 = produit.liste_sommets[j].id;
			if(clique[i] == 1 && clique[j] == 1 && depart[id1][id2] == 1 )
				somme++;
		}
	}
	//printf("%d valeur @@@@@@\n", somme);

	return somme;
}
void calcul_clique(int **matrice,int sommets,int *dans_clique,int taille_clique,int *candidat,int taille_candidat,GRAPHE_CYCLE produit, int **depart)
{	//calcul de la clique max recursif
	

	int i,j;
	
	if( taille_candidat == 0)
	{
		int val1, val2;
		//printf("herer %d. %d \n",taille_clique,taille_clique_max);
		val1 = somme_sommets_aretes(dans_clique_max,produit,depart, sommets);
		val2 = somme_sommets_aretes(dans_clique,produit,depart, sommets);
		if(val2 > val1)
		{
			taille_clique_max = taille_clique;
			for (i = 0 ;  i < sommets ; i++)
				dans_clique_max[i] = dans_clique[i];
		}
		return;
	}
	
	if (taille_candidat + taille_clique <= taille_clique_max)
	{
		//printf("je suis ici \n");
		return;
	}
	
	//else
	//if(taille_candidat == 1)printf("one time\n"); 
	int taille_candidat_temp;
	int *candidat_temp;
	candidat_temp = malloc( sommets * sizeof(int));


	/*if(taille_candidat == 1) {
		for (j = 0 ;  j < sommets ; j ++)
			{	printf("%d ", candidat[j]);}printf(" fin\n");
	}*/
	for (i = 0 ;  i < sommets ; i ++)
	{
		if ( candidat[i] == 1)
		{
			//if(taille_candidat == 1) printf("is %d \n",i);
			candidat[i] = 0;
			dans_clique[i] = 1 ;
			taille_candidat_temp = taille_candidat;
			
			for (j = 0 ;  j < sommets ; j ++)
			{
				candidat_temp[j] = candidat[j]; 
				if ((candidat[j] == 1) && (matrice[i][j] == 0))
				{
					candidat_temp[j] = 0;
					taille_candidat_temp--;	
				}	
			}
			
			taille_candidat_temp--;
			calcul_clique(matrice,sommets,dans_clique,taille_clique + 1,candidat_temp,taille_candidat_temp,produit,depart);
			dans_clique[i] = 0;
			candidat[i] = 1;
		}	
		
	}
	free(candidat_temp);
		 
}
void la_clique_max( int **matrice, int sommets,int **depart,GRAPHE_CYCLE produit)
{	//Debut calcul de la clique -- Initialisation
	int i;
	int *candidat;
	int *dans_clique;
	
	dans_clique_max = malloc( sommets *sizeof(int));
	if (!dans_clique_max) { fprintf(stderr,"cannot malloc dans_clique_max %d\n",sommets); exit(41); }
	dans_clique = malloc( sommets *sizeof(int));
	if (!dans_clique) { fprintf(stderr,"cannot malloc dans_clique %d\n",sommets); exit(42); }
	candidat = malloc( sommets *sizeof(int));
	if (!candidat) { fprintf(stderr,"cannot malloc candidat %d\n",sommets); exit(43); }
	
	//initialisation 
	for(i = 0; i < sommets ; i++ )
	{
		candidat[i] 	= 1;
		dans_clique[i]	= 0;
		dans_clique_max[i] = 0;
	}
	
	taille_clique_max = 0;

	
	calcul_clique(matrice,sommets,dans_clique,0,candidat,sommets,produit,depart); // 0 taille de la clique initial  et m.nb_atome = nb sommets candidats
	 
	free(dans_clique);
	free(candidat);
	
}

GRAPHE_CYCLE construction_commun_max(GRAPHE_CYCLE a,GRAPHE_CYCLE p)
{

	GRAPHE_CYCLE c;
	int tab[a.nb_sommets];
	int i, nb_at= 0;
	for(i=0;i < a.nb_sommets ; i++)
	{
		tab[i] = 0;
	}	

	for(i=0;i < p.nb_sommets ; i++)
	{
		
		if(dans_clique_max[i] == 1)
		{
			tab[p.liste_sommets[i].id] = 1;
		}
			
	}
	//free(dans_clique_max);

	for(i=0;i < a.nb_sommets ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}
	//printf("na atmes = %d\n",nb_at);
	int j,nb_liaisons = 0;
	for(i=0;i < a.nb_sommets - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < a.nb_sommets; j++)
			{
				if(tab[j] == 1 && poids_arete_graphe_cycle(a,a.liste_sommets[i].id,a.liste_sommets[j].id) != -1)
						nb_liaisons ++;
			}
		}
	}
	
	c.nb_aretes = nb_liaisons;
	c.nb_sommets = nb_at;
	//printf("na liaisons s = %d\n",nb_liaisons);
	return c;



}

void liberation_matrice(int **matrice, int taille)
{
	int i;
	for ( i = 0; i < taille ;i++)
		free(matrice[i]);
	free(matrice);
}

struct couple *construction_couples_cycles(GRAPHE_CYCLE a,GRAPHE_CYCLE b,int taille)
{//Construction des couples de cycles  compatibles
	
	int n = 0,i,j;
	
	struct couple *couple_at;
	couple_at = malloc(taille * sizeof(couple));

	for(i= 0; i < a.nb_sommets; i++) 
	{ 
		for(j= 0; j < b.nb_sommets;j++)
		{
			if( (float)(valeur_absolue(a.liste_sommets[i].poids - b.liste_sommets[j].poids)) <= 0.2 *min(a.liste_sommets[i].poids ,b.liste_sommets[j].poids))
			{
				couple_at[n].a1 = i;
				couple_at[n].a2 = j;
				n++;
			}
		}
	}
	//printf("%d %d", n , taille);
	return couple_at;

}

struct type_arete ** construction_matrice_graphe_cycles(GRAPHE_CYCLE a)
{//construction de la matrice de liaison d'une molecule
	
	int i,j;
	struct type_arete **matrice = NULL;

	matrice =  malloc(a.nb_sommets * sizeof(struct type_arete *));
	
	for(i=0; i < a.nb_sommets ;i++) matrice[i] =  malloc(a.nb_sommets * sizeof(struct type_arete));

	if(matrice  == NULL)
		probleme_memoire();


	for(i=0; i < a.nb_sommets ;i++)
	{
		for(j=0;j< a.nb_sommets;j++)
		{	
			matrice[i][j].type = AUCUNE_LIAISON;
			matrice[i][j].poids = AUCUNE_LIAISON;
			//printf("(%d,%d)", matrice[i][j].type,matrice[i][j].poids);
		}
		//printf("\n");
	}
	//printf("\n after \n");
	for(i =0; i< a.nb_aretes;i++)
	{
		matrice[a.liste_aretes[i].id1][a.liste_aretes[i].id2].type= a.liste_aretes[i].type;
		matrice[a.liste_aretes[i].id1][a.liste_aretes[i].id2].poids= a.liste_aretes[i].poids;
		matrice[a.liste_aretes[i].id2][a.liste_aretes[i].id1].type= a.liste_aretes[i].type;
		matrice[a.liste_aretes[i].id2][a.liste_aretes[i].id1].poids= a.liste_aretes[i].poids;
	}

	/*for(i=0; i < a.nb_sommets ;i++)
	{
		for(j=0;j< a.nb_sommets;j++)
		{	
			printf("(%d,%d)", matrice[i][j].type,matrice[i][j].poids);
		}
		printf("\n");
	}*/


	return matrice;
	
}

void liberer_type_arete( struct type_arete **m, GRAPHE_CYCLE a)
{
	int i;
	for(i= 0; i < a.nb_sommets; i++) free(m[i]);
	free(m);

	
}
graph graphe_produit_cycles(GRAPHE_CYCLE a,GRAPHE_CYCLE b)
{ //prend en entrée deux graphes de cycles  et contruit le graphe produit
		
		
	
	int taille = 0;
	
	//calcul de la taille du graphe produit
	int i,j;
	for(i= 0; i < a.nb_sommets; i++) 
		for(j= 0; j < b.nb_sommets;j++)
			if( (float)(valeur_absolue(a.liste_sommets[i].poids - b.liste_sommets[j].poids)) <= 0.2 *min(a.liste_sommets[i].poids ,b.liste_sommets[j].poids))
				taille ++;
	
	//couple de liaisons entre les nouveaux sommets


	struct couple * couple_atome = construction_couples_cycles(a,b,taille);

	
	//construction de la matrice de liaison d'une molecule
	struct type_arete **m1 = construction_matrice_graphe_cycles(a);
	struct type_arete **m2 =construction_matrice_graphe_cycles(b);

   
  	int** matrice_liaisons = NULL;

	if( taille > 0)
	{	
		matrice_liaisons =  malloc(taille * sizeof(int *));	
		for(i = 0 ;i < taille;i++) 
    		matrice_liaisons[i] =  malloc(taille * sizeof(int));
	}
	
  //remplissage des liaisons
	int i1,i2,j1,j2;
	for(i= 0; i < taille ;i++)
		for(j= 0; j < taille ;j++)
			matrice_liaisons[i][j] = 0;

	for(i= 0; i < taille ;i++)
		matrice_liaisons[i][i] = 1;

	for(i= 0; i < taille ;i++){ 
		i1 = couple_atome[i].a1;
		i2 = couple_atome[i].a2;
		for(j= i + 1 ; j < taille ;j++){
			j1=couple_atome[j].a1;
			j2=couple_atome[j].a2;
			
			if( m1[i1][j1].type == m2[i2][j2].type && m1[i1][j1].poids == m2[i2][j2].poids){
				if(((i1 == j1) && (m2[i2][j2].type != AUCUNE_LIAISON)) || ((i2 == j2) && (m1[i1][j1].type != AUCUNE_LIAISON)) ||( (i1 != j1) && (i2!=j2))||( (i1 ==j1) && (i2==j2)) )
				{
					matrice_liaisons[i][j] = 1;
					matrice_liaisons[j][i] = 1;					
				}
			}
		}
	}
	if(couple_atome != NULL)
		free(couple_atome);
	liberer_type_arete(m1,a);
	liberer_type_arete(m2,b);
	
	graph f = build_graph_from_matrix(taille, matrice_liaisons);	

	return f;
}

int*  graphe_g12_cycles(graph g12, int* clique_max,GRAPHE_CYCLE a,GRAPHE_CYCLE b)
{ //contruction du graphe commun

	int* taille_graphe_commun = (int*) malloc(sizeof(int)*2);
  	int nb_at =0,nb_liaisons=0,i,j,i1,j1;
	
  	int taille = nbnodes(g12);
	struct couple * couple_atome = construction_couples_cycles(a,b,taille);
	struct type_arete **m1 = construction_matrice_graphe_cycles(a);
	
	int tab[taille];
	for(i=0;i < taille ; i++)
	{
		tab[i] = clique_max[i];
	}
	
	for(i=0;i < taille - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j=i+1;j < taille ; j++)
			{
				if(tab[j]==1 && (couple_atome[i].a1 == couple_atome[j].a1))
						tab[j] = 0;
			}
		}
	}

	for(i=0;i < taille ; i++)
	{
		if(tab[i] ==1)
			nb_at++;
	}

	for(i=0;i < taille - 1; i++)
	{
		if(tab[i] == 1)
		{
			for(j = i+1;j < taille; j++)
			{
				if(tab[j] == 1)
				{
					
					i1 = couple_atome[i].a1;
					j1 = couple_atome[j].a1;
					if(m1[i1][j1].type != AUCUNE_LIAISON)
						nb_liaisons ++;
				}
			}
		}
	}
	free(couple_atome);
	liberer_type_arete(m1,a);
	taille_graphe_commun[0] = nb_at;
	taille_graphe_commun[1] = nb_liaisons;
	return taille_graphe_commun;
	
}

float similarite(GRAPHE_CYCLE a,GRAPHE_CYCLE b)
{

	
	//int 
	float sim = 0.0;
	if(a.nb_sommets == 0 || b.nb_sommets == 0)
		sim = -1.0;
	else
	{
		//comence ici
		graph g = graphe_produit_cycles(a,b);

		int taille = nbnodes(g);
		int** liaisons = build_matrix_from_graph(g);
		if(taille == 0)
		  	return 0;

		int* clique = clique_max(g, (long)date);

		int* taille_graphe_commun = graphe_g12_cycles(g,clique,a,b);

		int nb_atomes_communs = taille_graphe_commun[0];	
	  	int nb_laisons_communs = taille_graphe_commun[1];


	  	float num = (float)((nb_atomes_communs + nb_laisons_communs)*(nb_atomes_communs + nb_laisons_communs));
	  	float denum = (float)((a.nb_sommets + a.nb_aretes)*(b.nb_aretes + b.nb_sommets));
	  	sim = num/denum;
		free(taille_graphe_commun);
		free(clique);
		//fin ici
		
		if(dans_clique_max != NULL)
			free(dans_clique_max);
		taille_clique_max = 0;

		destroy(g);
	}	
	return sim;
}
