#include "motor.h"

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

void motorComputeMagneticPotential(motor *theMotor)
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
   
    femFullSystem* FullSystemStator = femFullSystemCreate(nElemDomain[0]); // créer un système de taille du nombre de noeud pour le maillage du stator (de numéro 0)



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
    return;
}


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



//
// ========= Projet à réaliser ===================
//