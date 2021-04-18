#include "motor.h"


/** motorMeshLocal
 * @in 
 * - const motorMesh* theMesh = ensemble des maillages du moteur considéré
 * - const int i = indice du ième sous-triangle considéré
 * - int* map = pointeur vers un tableau de 3 éléments préalablement alloué
 * - double* x = pointeur vers un tableau de 3 éléments préalablement alloué
 * - double* y = pointeur vers un tableau de 3 éléments préalablement alloué
 * @out 
 * rempli:
 * -> le tableau map par les 3 numéros des noeuds constituant le ième sous-triangle
 * -> le tableau x par les 3 abscisses des noeuds constituant le ième sous-triangle
 * -> le tableau y par les 3 ordonnées des noeuds constituant le ième sous-triangle
 */
void motorMeshLocal(const motorMesh* theMesh, const int iElem, int *map, double *x, double *y)
{    
    for (int j=0; j < 3; ++j) 
    {
        map[j] = theMesh->elem[iElem*3+j];
        x[j]   = theMesh->X[map[j]];
        y[j]   = theMesh->Y[map[j]]; 
    }   
}



/** renumberNodesDomain
 * @in 
 * - const motorMesh* theMesh = ensemble des maillages du moteur considéré
 * @out 
 * - pointeur vers un tableau contenant les noeuds triés par ordre d'appartenance à leur sous-domaine
 *  si les domaine 0, 1, ..., 11 contiennent respectivement 1, 2, ... 12 triangles chacun,
 *  alors le tableau contiendra d'abord les 3 noeuds du seul triangle du domaine 1, puis les 6 (2*3) noeuds du domaine 2
 *  et ainsi de suite jusqu'au domaine 11 et ses 36 (12*3) noeuds. 
 */
int* renumberNodesDomain(const motorMesh* theMesh)
{

    int totalNodes = 0;
    for(int i =0 ; i< theMesh->nDomain; i++)
    {
        totalNodes += theMesh->nElemDomain[i];
    }
    totalNodes = totalNodes * 3;

    int* mapping = malloc(sizeof(int) * totalNodes);
    
    int indexDomain0 = 0;
    int indexDomain1 = indexDomain0 + theMesh->nElemDomain[0]*3;
    int indexDomain2 = indexDomain1 + theMesh->nElemDomain[1]*3;
    int indexDomain3 = indexDomain2 + theMesh->nElemDomain[2]*3;
    int indexDomain4 = indexDomain3 + theMesh->nElemDomain[3]*3;
    int indexDomain5 = indexDomain4 + theMesh->nElemDomain[4]*3;
    int indexDomain6 = indexDomain5 + theMesh->nElemDomain[5]*3;
    int indexDomain7 = indexDomain6 + theMesh->nElemDomain[6]*3;
    int indexDomain8 = indexDomain7 + theMesh->nElemDomain[7]*3;
    int indexDomain9 = indexDomain8 + theMesh->nElemDomain[8]*3;
    int indexDomain10 = indexDomain9 + theMesh->nElemDomain[9]*3;
    int indexDomain11 = indexDomain10 + theMesh->nElemDomain[10]*3;

    //printf("indexDomain0 : %d \n",indexDomain0);
    //printf("indexDomain1 : %d \n",indexDomain1);
    //printf("indexDomain2 : %d \n",indexDomain2);
    //printf("indexDomain3 : %d \n",indexDomain3);
    //printf("indexDomain4 : %d \n",indexDomain4);
    //printf("indexDomain5 : %d \n",indexDomain5);
    //printf("indexDomain6 : %d \n",indexDomain6);
    //printf("indexDomain7 : %d \n",indexDomain7);
    //printf("indexDomain8 : %d \n",indexDomain8);
    //printf("indexDomain9 : %d \n",indexDomain9);
    //printf("indexDomain10 :%d  \n",indexDomain10);
    //printf("indexDomain11 :%d  \n",indexDomain11);

    for(int iTriangle =0; iTriangle < theMesh->nElem; iTriangle++)
    //Parcourir tous les triangles
    {
        if( theMesh->domain[iTriangle] == 0 ) // correspond au domaine du stator
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain0] = theMesh->elem[iTriangle*3];
            mapping[indexDomain0+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain0+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain0],mapping[indexDomain0+1],mapping[indexDomain0+2]); 
            indexDomain0 +=3;
        }
        else if( theMesh->domain[iTriangle] == 1 ) // correspond au domaine COIL_AP
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain1] = theMesh->elem[iTriangle*3];
            mapping[indexDomain1+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain1+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain1],mapping[indexDomain1+1],mapping[indexDomain1+2]); 
            indexDomain1 +=3;
        }
        else if( theMesh->domain[iTriangle] == 2 ) // correspond au domaine COIL_AN
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain2] = theMesh->elem[iTriangle*3];
            mapping[indexDomain2+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain2+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain2],mapping[indexDomain2+1],mapping[indexDomain2+2]); 
            indexDomain2 +=3;
        }
        else if( theMesh->domain[iTriangle] == 3 ) // correspond au domaine COIL_BP
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain3] = theMesh->elem[iTriangle*3];
            mapping[indexDomain3+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain3+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain3],mapping[indexDomain3+1],mapping[indexDomain3+2]); 
            indexDomain3 +=3;
        }
        else if( theMesh->domain[iTriangle] == 4 ) // correspond au domaine COIL_BN
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain4] = theMesh->elem[iTriangle*3];
            mapping[indexDomain4+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain4+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain4],mapping[indexDomain4+1],mapping[indexDomain4+2]); 
            indexDomain4 +=3;
        }
        else if( theMesh->domain[iTriangle] == 5 ) // correspond au domaine COIL_CP
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain5] = theMesh->elem[iTriangle*3];
            mapping[indexDomain5+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain5+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain5],mapping[indexDomain5+1],mapping[indexDomain5+2]); 
            indexDomain5 +=3;
        }
        else if( theMesh->domain[iTriangle] == 6 ) // correspond au domaine COIL_CN
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain6] = theMesh->elem[iTriangle*3];
            mapping[indexDomain6+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain6+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain6],mapping[indexDomain6+1],mapping[indexDomain6+2]); 
            indexDomain6 +=3;
        }
        else if( theMesh->domain[iTriangle] == 7 ) // correspond au domaine Stator_gap
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain7] = theMesh->elem[iTriangle*3];
            mapping[indexDomain7+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain7+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain7],mapping[indexDomain7+1],mapping[indexDomain7+2]); 
            indexDomain7 +=3;
        }
        else if( theMesh->domain[iTriangle] == 8 ) // correspond au domaine Rotor_core
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain8] = theMesh->elem[iTriangle*3];
            mapping[indexDomain8+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain8+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain8],mapping[indexDomain8+1],mapping[indexDomain8+2]); 
            indexDomain8 +=3;
        }
        else if( theMesh->domain[iTriangle] == 9 ) // correspond au domaine Rotor_air
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain9] = theMesh->elem[iTriangle*3];
            mapping[indexDomain9+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain9+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain9],mapping[indexDomain9+1],mapping[indexDomain9+2]); 
            indexDomain9 +=3;
        }
        else if( theMesh->domain[iTriangle] == 10 ) // correspond au domaine Rotor_gap
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain10] = theMesh->elem[iTriangle*3];
            mapping[indexDomain10+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain10+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain10],mapping[indexDomain1O+1],mapping[indexDomain10+2]); 
            indexDomain10 +=3;
        }
        else if( theMesh->domain[iTriangle] == 11 ) // correspond au domaine Air_gap
        {
            // recopier les 3 noeuds dans le mapping;
            mapping[indexDomain11] = theMesh->elem[iTriangle*3];
            mapping[indexDomain11+1] = theMesh->elem[iTriangle*3+1];
            mapping[indexDomain11+2] = theMesh->elem[iTriangle*3+2];
            //printf("%d : %d %d %d\n",iTriangle, mapping[indexDomain11],mapping[indexDomain11+1],mapping[indexDomain11+2]); 
            indexDomain11 +=3;
        }
        else
        {
            printf("An error occured in renumbering nodes \n"); 
        }

        //printf("%d : %d %d %d\n",iTriangle, theMesh->elem[iTriangle*3],theMesh->elem[iTriangle*3+1],theMesh->elem[iTriangle*3+2]); 


    }
    printf("No problem until here \n");

        	

    return mapping;
}

void* createStatorMesh(const motorMesh* theMesh, const femMesh* statorMesh)
{
    
}
//
// ========= Projet à réaliser ===================
//

void motorAdaptMesh(motor *theMotor, double delta)
{
    motorMesh *theMesh = theMotor->mesh;
    
    double x,y;
    for(int i = 0; i < theMesh->nNode; ++i){
        if  (theMotor->movingNodes[i] == 1){
            x = theMesh->X[i]*cos(delta) - theMesh->Y[i]*sin(delta);
            y = theMesh->X[i]*sin(delta) + theMesh->Y[i]*cos(delta);
            theMesh->X[i] = x;
            theMesh->Y[i] = y; }}
    theMotor->theta += delta;
}

double motorComputeCouple(motor *theMotor)
{
    return 0.0;

}

void motorComputeCurrent(motor *theMotor)
{
    return;
}

void motorComputeMagneticPotential(motor* theMotor)
{
    // Résolution inspirée du devoir 4 ainsi que de sa correction disponible sur https://www.youtube.com/watch?v=580gEIVVKe8.
   
    int size = theMotor->size; 
    motorMesh* theMesh = theMotor->mesh;

    double* a = theMotor->a;
    double* js = theMotor->js;
    double* mu = theMotor->mu;

    int* elem = theMesh->elem;
    double* X = theMesh->X;
    double* Y = theMesh->Y;
    int nElem = theMesh->nElem;
    int nNode = theMesh->nNode;
    int nLocalNode = theMesh->nLocalNode; // = 3 car tjrs un triangle 
    int nDomain = theMesh->nDomain;
    int* nElemDomain = theMesh->nElemDomain;
    string* nameDomain = theMesh->nameDomain;
    int* domain = theMesh->domain;

    /*
    femMesh* statorMesh = malloc(sizeof(femMesh));
    statorMesh->elem = elem[0]; // le stator est le premier
    statorMesh->X = &X[0];
    statorMesh->Y = &Y[0];
    statorMesh->nElem = &nElem[0];
    statorMesh->nNode = &nNode[0];
    statorMesh->nLocalNode = 3;
    //statorMesh->number = ??;
    */
   
    int nNodeStator = nElemDomain[0]*3;
    femFullSystem* FullSystemStator = femFullSystemCreate(nNodeStator); // créer un système de taille du nombre de noeud pour le maillage du stator (de numéro 0)
    printf("nNodeStator : %d \n", nNodeStator);

    for(int iTriangle =0; iTriangle< nElem; iTriangle++)
    // Parcourir tous les triangles des maillages
    {
        if( strcmp( nameDomain[domain[iTriangle]], "Rotor_core") == 0)
        // Imposer un potentiel magnétique constant et nul dans le rotor
        {   
            //printf(" Le triangle %d est formé des 3 noeuds : %d %d %d \n", iTriangle, elem[iTriangle*3], elem[iTriangle*3+1], elem[iTriangle*3+2]);
            //printf(" le triangle %d appartient au sous-domaine: %d intitulé:  %s\n", iTriangle, domain[iTriangle], nameDomain[domain[iTriangle]] );
            a[elem[iTriangle*3]] = 0.0;
            a[elem[iTriangle*3+1]] = 0.0;
            a[elem[iTriangle*3+2]] = 0.0;
        }
        else if(strcmp( nameDomain[domain[iTriangle]], "Stator_core") == 0) 
        // Résoudre l'équation de Poisson dans le stator
        {
            // 1. Construire le système linéaire = assembler les matrices locales

            // 1.1. Récupérer les numéros et coordonnées des noeuds du sous-triangle courant 
            int* map = malloc(sizeof(int) * nLocalNode); if(map == NULL) {Error("Error malloc map in motorComputeMagneticPotential");}
            double* x = malloc(sizeof(double) * nLocalNode); if(map == NULL) {Error("Error malloc x in motorComputeMagneticPotential");}
            double* y = malloc(sizeof(double) * nLocalNode); if(map == NULL) {Error("Error malloc y in motorComputeMagneticPotential");}
            motorMeshLocal(theMesh, iTriangle, map, x, y); 

            // 1.2 Calculer l'intégrale du sous-triangle avec la règle d'intégration 
            // de Hammer à un point
            
            // Effectuer l'intégrale sur le repère du parent
            double xsi = 0.333333333333333;
            double eta = 0.333333333333333;
            double weight = 1.0;

            // fonctions de forme dans le triangle parent
            double* phi = malloc(sizeof(double) * 3); // car 3 sommets à chaque triangle
            phi[0] = 1.0 - xsi - eta;   
            phi[1] = xsi;
            phi[2] = eta;      
            //#femDiscretePhi2       
            
            // dérivées des fonctions de formes dans le triangle parent
            double* dphidxsi = malloc(sizeof(double)*3);
            double* dphideta = malloc(sizeof(double)*3);
            dphidxsi[0] = -1.0;  
            dphidxsi[1] =  1.0;
            dphidxsi[2] =  0.0;
            dphideta[0] = -1.0;  
            dphideta[1] =  0.0;
            dphideta[2] =  1.0;
            //#femDiscreteDphi2(space,xsi,eta,dphidxsi,dphideta); 
            

            // Interpolation par des fonctions de forme
            // calcul du gradient de la transformation
            double dxdxsi = 0.0;                    double dydxsi = 0.0; 
            double dxdeta = 0.0;                    double dydeta = 0.0;
            
            for(int i = 0; i < 3 ; i++) 
            {
                dxdxsi += x[i] * dphidxsi[i];           dydxsi += y[i] * dphidxsi[i];     
                dxdeta += x[i] * dphideta[i];           dydeta += y[i] * dphideta[i]; 
            }
            
            
            // calcul du Jacobien de la transformation
            double Jacobien = dxdxsi * dydeta - dxdeta * dydxsi;
            if(Jacobien <0) // noeuds sont mal orientés => nécessité de réorienter les noeuds dans le maillage
            {
                int node = elem[nLocalNode*iTriangle];
                elem[nLocalNode*iTriangle] = elem[nLocalNode*iTriangle+2];
                elem[nLocalNode*iTriangle + 2] = node;
            }
            Jacobien = fabs(Jacobien);

            // calcul du gradient des fonctions de forme par interpolation des fonctions de formes
            double* dphidx = malloc(sizeof(double) * 3);
            double* dphidy = malloc(sizeof(double) * 3);
            for (int i = 0; i < 3; i++) 
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / Jacobien ;      
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / Jacobien ; 
            }       

            // Assemblage des matrices locales dans la matrice globale
            for (int i = 0; i < 3; i++) 
            { 
                FullSystemStator->B[map[i]] += weight * Jacobien * phi[i];// Assemblage du vecteur B
                for( int j = 0; j < 3; j++) // Assemblage de la matrice A
                {
                    /*
                    printf("Strange values ...\n");
                    printf("%d \n", FullSystemStator->A[map[i]][map[j]]);
                    */
                    FullSystemStator->A[map[i]][map[j]] += weight * Jacobien * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]);
                }
            }
    
            

            free(map);
            free(x);
            free(y);

            free(phi);
            free(dphidxsi);
            free(dphideta);

            free(dphidx);
            free(dphidy);

            
        }
    }


    // ARRIVER ICI : IMPOSER LA CONTRAINTE SUR LE RAYON EXTERIEUR DU STATOR

    // Imposer la condition limite pour le stator: 
    // potentiel nul sur son rayon extérieur

    // Déterminer la valeur du rayon extérieur
    double rayonExt = 0.0;
    double rayon = 0.0;
    int nb = 0;
    for(int i=0; i< theMesh->nElemDomain[0]; i++)
    //Parcourir chaque noeud du stator et calculer son rayon
    {
        rayon = sqrt( X[i]*X[i] + Y[i]*Y[i] );
        //printf("Rayon : %f \n", rayon);

        if( rayon >= rayonExt)
        {
            nb++;
            rayonExt = rayon;
        }
    }
    printf("Rayon extérieur : %f \n", rayonExt);
    printf("Nombre de noeuds sur le rayon extérieur : %d \n", nb);

    // Imposer une valeur nulle pour les noeuds situés sur le rayon extérieur.
    for(int i=0; i< theMesh->nElemDomain[0]; i++)
    //Parcourir chaque noeud du stator et calculer son rayon
    {
        rayon = sqrt( X[i]*X[i] + Y[i]*Y[i] );
        if( rayon == rayonExt)
        {
            a[elem[i]] = 0.0; 
        }
    }

    /*
    femEdges* edges = femEdgesCreate(theMesh);
    for (int iEdge = 0; iEdge < edges->nEdge; iEdge++) 
    {      
        femEdge* edge = &(edges->edges[iEdge]);
        if (edge->elem[1] == -1) // détecter les frontières pour y imposer la condition limite
        {  
            //imposer la contrainte
            double value = 0.0; // valeur de la condition limite
            for (int j = 0; j < 2; j++) 
            {
                femFullSystemConstrain(theProblem->system, edge->node[j] ,value);  
            }
        }
            
    }
    
    femFullSystemEliminate(FullSystemStator);
    // modifier la valeur dans le tableau a
    
    // Ne pas oublier tous les free !!!
    femFullSystemFree(FullSystemStator);
    */
    
    return;

} 

//
// ========= Projet à réaliser ===================
//