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
    int nDomain = theMesh->nDomain;
    int* nElemDomain = theMesh->nElemDomain;
    string* nameDomain = theMesh->nameDomain;
    int* domain = theMesh->domain;

    
    for(int iElem =0; iElem< nElem; iElem++)
    {
        if( strcmp( nameDomain[domain[iElem]], "Rotor_core") == 0)
        // Imposer un potentiel magnétique constant et nul dans le rotor
        {   
            //printf(" Le triangle %d est formé des 3 noeuds : %d %d %d \n", iElem, elem[iElem*3], elem[iElem*3+1], elem[iElem*3+2]);
            //printf(" le triangle %d appartient au sous-domaine: %d intitulé:  %s\n", iElem, domain[iElem], nameDomain[domain[iElem]] );
            a[elem[iElem*3]] = 0.0;
            a[elem[iElem*3+1]] = 0.0;
            a[elem[iElem*3+2]] = 0.0;
        }
    }

    return;
}

/*

        else if(strcmp( nameDomain[domain[iElem]], "Stator_core") == 0) 
        // Résoudre l'équation de Poisson dans le stator
        {
            // Récupérer la liste des triangles dans le maillage stator_core

            // utiliser femMeshLocal


        }
*/

//
// ========= Projet à réaliser ===================
//
/** femMeshLocal
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
void femMeshLocal(const motorMesh* theMesh, const int i, int *map, double *x, double *y)
{
    for(int j = 0; j < theMesh->nLocalNode; j++)
    {
        map[j] = theMesh->elem[i * theMesh->nLocalNode + j];
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
    }
}


void femPoissonSolve(femPoissonProblem *theProblem)
{
    femMesh* mesh = theProblem->mesh;
    femIntegration* rule = theProblem->rule;
    femDiscrete * space = theProblem->space;

    int nLocalNode = mesh->nLocalNode; // un triangle ou un quadrilatère 

    double** A = theProblem->system->A;
    double* B = theProblem->system->B;

    // 1. Construire le système linéaire = assembler les matrices locales

   
    for(int iElem=0; iElem< theProblem->mesh->nElem; iElem++)  // Pour chaque élément local (sous-triangle ou sous-quadrilatère) du maillage
    {
        // 1.1. Récupérer les indices et valeurs des noeuds (via la fonction femMeshLocal précédente)
                
        // coordonnées des noeuds du sous-triangle ou sous-quadrilatère
        double* x = malloc(sizeof(double) * mesh->nLocalNode);      
        double* y = malloc(sizeof(double) * mesh->nLocalNode);
        // indices des noeuds du sous-triangle ou sous-quadrilatère
        int* map = malloc(sizeof(int) * mesh->nLocalNode);

        
        femMeshLocal(theProblem -> mesh, iElem, map, x, y); //map, x, y ont été remplis


        // 1.2 Calculer l'intégrale du sous-triangle/ sous-quadrilatère avec la
        //     règle d'intégration fournie

        for(int iIntegral = 0; iIntegral < rule->n; iIntegral++)
        {
            // Effectuer l'intégrale sur le repère du parent
            double xsi = rule->xsi[iIntegral];
            double eta = rule->eta[iIntegral];
            double weight = rule->weight[iIntegral];

            // fonctions de forme dans le triangle parent
            double* phi = malloc(sizeof(double) * space->n);
            femDiscretePhi2(space,xsi, eta, phi); // phi a été remplie

            // dérivées des fonctions de formes dans le triangle parent
            double* dphidxsi = malloc(sizeof(double)*space->n);
            double* dphideta = malloc(sizeof(double)*space->n);
            femDiscreteDphi2(space,xsi,eta,dphidxsi,dphideta); // dphidxsi et dphideta ont été remplies


            // calcul du gradient de la transformation
            double dxdxsi = 0.0;                    double dydxsi = 0.0; 
            double dxdeta = 0.0;                    double dydeta = 0.0;
            
            for(int i = 0; i < space->n ; i++)
            {
                dxdxsi += x[i] * dphidxsi[i];           dydxsi += y[i] * dphidxsi[i];     
                dxdeta += x[i] * dphideta[i];           dydeta += y[i] * dphideta[i]; 
            }
            
            // calcul du Jacobien de la transformation
            double Jacobien = dxdxsi * dydeta - dxdeta * dydxsi;
            if(Jacobien <0) // noeuds sont mal orientés => nécessité de réorienter les noeuds dans le maillage
            {
                int node = mesh->elem[nLocalNode*iElem];
                mesh->elem[nLocalNode*iElem] = mesh->elem[nLocalNode*iElem+2];
                mesh->elem[nLocalNode*iElem + 2] = node;
            }
            Jacobien = fabs(Jacobien);

            // calcul du gradient des fonctions de forme par interpolation des fonctions de formes
            double* dphidx = malloc(sizeof(double) * space->n);
            double* dphidy = malloc(sizeof(double) * space->n);
            for (int i = 0; i < space->n; i++) 
            {    
                dphidx[i] = (1 / Jacobien) * (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) ;      
                dphidy[i] = (1 / Jacobien) * (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) ; 
            }       

            // Assemblage des matrices locales dans la matrice globale
            for (int i = 0; i < space->n; i++) 
            { 
                B[map[i]] += (phi[i] * Jacobien) * weight;// Assemblage du vecteur B
                
                for( int j = 0; j < space->n; j++) // Assemblage de la matrice A
                {
                    A[map[i]][map[j]] += dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j] * Jacobien * weight;
                }
            }
            
            free(dphidx);
            free(dphidy);
            free(dphidxsi);
            free(dphideta);
            free(phi);
        }

        free(x);
        free(y);
        free(map);
        
    }     

    // Imposer la contrainte
    femEdges* edges = theProblem->edges;
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

    femFullSystemEliminate(theProblem -> system);
}
