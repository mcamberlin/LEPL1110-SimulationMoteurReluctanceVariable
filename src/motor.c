#include "motor.h"
#include <float.h>

static motorMesh* theGlobalMotorMesh;

typedef enum 
{
    Stator_core = 0,
    Rotor_core = 8,   
    Coil_AP = 1,
    Coil_AN = 2,
    Coil_BP = 3,
    Coil_BN = 4,
    Coil_CP = 5,
    Coil_CN = 6,
    Rotor_gap = 10,
    Air_gap = 11,
    Stator_gap = 7,
    Rotor_air = 9
} iDomain;

void printFemMesh(const femMesh* theFemMesh)
{
    int start = 0;
    int stop = 10;
    printf("\n ====== femMesh information ============================\n");
    printf("    elem = [ ");
    for(int i = start; i < stop; i++)
    {
        printf("%d, ",theFemMesh->elem[i]); 
    }
    printf("... ]\n");

    printf("    X = [ ");
    for(int i = start ;i < stop; i++)
    {
        printf("%f, ",theFemMesh->X[i]);
    }
    printf("... ]\n");

    printf("    Y = [ ");
    for(int i = start ;i < stop; i++)
    {
        printf("%f, ",theFemMesh->Y[i]);
    }
    printf("... ]\n");
    printf("    Number of elements           : %d\n",theFemMesh->nElem);   
    printf("    Number of nodes              : %d\n",theFemMesh->nNode);  
    printf("    Number of local nodes        : %d\n",theFemMesh->nLocalNode); 
        
    printf("    number = [ ");
    for(int i = start ;i < stop; i++)
    {
        printf("%d, ",theFemMesh->number[i]);
    }
    printf("... ]\n");
    printf("=========================================================================\n");

}

void printTriangle(const motorMesh* theMotorMesh, const int iTriangle)
{
    double* X = theMotorMesh->X;
    double* Y = theMotorMesh->Y;
    int* elem = theMotorMesh->elem;

    printf("======= Triangle %d =======\n", iTriangle);    
    for(int i =0; i< 3; i++)
    {
        printf("\t  %d : ( %le, %le )\n", elem[iTriangle*3+i], X[elem[iTriangle*3+i]], Y[elem[iTriangle*3+i]]);
    }
    printf("\n");    
}   

void printDoubleArray(const char* name, const double* ptr)
{
    int start = 0;
    int stop = 10;
    printf("\n ==================================\n");
    printf(" %s = [ ", name);
    for(int i = start; i < stop; i++)
    {
        printf("%f, ",ptr[i]); 
    }
    printf("... ]\n");
}

void printIntArray(const char* name, const int* ptr)
{
    int start = 0;
    int stop = 10;
    printf("\n ==================================\n");
    printf(" %s = [ ", name);
    for(int i = start; i < stop; i++)
    {
        printf("%d, ",ptr[i]); 
    }
    printf("... ]\n");
}
//    
// ========= Fonctions utiles motorComputeMagneticPotential() ===================
//
    /** motorMeshToFemMeshConverter
     * @in 
     * - const motorMesh* theMesh = ensemble des maillages du moteur considéré
     * @out 
     * - femMesh* = structure femMesh* contenant les mêmes attributs que ceux dans motorMesh 
     * (à l'exception de ceux liés aux sous-domaines). 
     */
    femMesh* motorMeshToFemMeshConverter(const motorMesh* theMotorMesh)
    {
        
        femMesh* theFemMesh = malloc(sizeof(femMesh));
        theFemMesh->elem = theMotorMesh->elem;
        theFemMesh->X = theMotorMesh->X;
        theFemMesh->Y = theMotorMesh->Y;
        theFemMesh->nElem = theMotorMesh->nElem;
        theFemMesh->nNode = theMotorMesh->nNode;
        theFemMesh->nLocalNode = theMotorMesh->nLocalNode;
        theFemMesh->number = malloc(sizeof(int)*theFemMesh->nNode); //Utilisé pour la renumérotation des noeuds dans le maillage 
        for (int i = 0; i < theFemMesh->nNode; i++)
        {
            theFemMesh->number[i] = i; 
        } 
        return theFemMesh;
    }

    void freeFemMeshConverted(femMesh* theFemMesh)
    {
        free(theFemMesh->number);
        free(theFemMesh);
    }

    /** femMeshLocal
     * @in 
     * - const femMesh* theMesh = maillage considéré
     * - const int iElem = indice du ième sous-triangle considéré
     * - int* map = pointeur vers un tableau de 3 éléments préalablement alloué
     * - double* x = pointeur vers un tableau de 3 éléments préalablement alloué
     * - double* y = pointeur vers un tableau de 3 éléments préalablement alloué
     * @out 
     * rempli:
     * -> le tableau map par les 3 numéros des noeuds constituant le ième sous-triangle
     * -> le tableau x par les 3 abscisses des noeuds constituant le ième sous-triangle
     * -> le tableau y par les 3 ordonnées des noeuds constituant le ième sous-triangle
     */
    void myFemMeshLocal(const femMesh *theMesh, const int iElem, int *map, double *x, double *y)
    {
        int j,nLocal = theMesh->nLocalNode;
        
        for (j=0; j < nLocal; ++j) 
        {
            map[j] = theMesh->elem[iElem*nLocal+j];
            x[j]   = theMesh->X[map[j]];
            y[j]   = theMesh->Y[map[j]]; 
        }   
    }

//
// ========= Fonctions utiles motorAdaptMesh() ===================
//

    int startAirGap()
    {
        int start = 0;
        for(int i=0; i< Air_gap; i++)
        {
            start += theGlobalMotorMesh->nElemDomain[i];
        }
        return start; 
    }

    int endAirGap()
    {
        return startAirGap() + theGlobalMotorMesh->nElemDomain[Air_gap];    
    }

    double radius(const int node)
    {
        double x = theGlobalMotorMesh->X[node];
        double y = theGlobalMotorMesh->Y[node];
        return sqrt(x*x + y*y);
    }

    double distanceBetweenNodes(const int node1, const int node2)
    {
        double x1 = theGlobalMotorMesh->X[node1];
        double y1 = theGlobalMotorMesh->Y[node1];
        double x2 = theGlobalMotorMesh->X[node2];
        double y2 = theGlobalMotorMesh->Y[node2];

        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    }

    double squarredDistanceBetweenNodes(const int node1, const int node2)
    {
        double x1 = theGlobalMotorMesh->X[node1];
        double y1 = theGlobalMotorMesh->Y[node1];
        double x2 = theGlobalMotorMesh->X[node2];
        double y2 = theGlobalMotorMesh->Y[node2];

        return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);   
    }



//
// ========= Projet à réaliser ===================
//

void motorAdaptMesh(motor *theMotor, double delta)
{
    motorMesh *theMotorMesh = theMotor->mesh;  
    theGlobalMotorMesh = theMotorMesh; 
    
    double* X = theMotorMesh->X;
    double* Y = theMotorMesh->Y;
    int* elem = theMotorMesh->elem;
    double L = theMotor->L;
    
    double x,y;
    for(int i = 0; i < theMotorMesh->nNode; ++i)
    {
        if  (theMotor->movingNodes[i] == 1)
        {
            x = theMotorMesh->X[i]*cos(delta) - theMotorMesh->Y[i]*sin(delta);
            y = theMotorMesh->X[i]*sin(delta) + theMotorMesh->Y[i]*cos(delta);   
            X[i] = x;
            Y[i] = y; 
        }
    }
    theMotor->theta += delta;
    
    
    // 1. Déterminer si les noeuds sur le rayon intérieur ou sur le rayon extérieur
        int startTriangleAirGap = startAirGap();
        int endTriangleAirGap = endAirGap();
        //printf("startAirGap : %d, endAirGap : %d \n", startTriangleAirGap, endTriangleAirGap);

        //printTriangle(theMotorMesh,startTriangleAirGap);

        // calculer de la largeur annulaire : calculer le rayon de 3 noeuds d'un triangle dans Air_gap
        
        double rayon1 = radius(elem[startTriangleAirGap*3]);
        double rayon2 = radius(elem[startTriangleAirGap*3+1]);
        double rayon3 = radius(elem[startTriangleAirGap*3+2]);
        double rayonMax = fmax(rayon1, fmax(rayon2, rayon3));
        double rayonMin = fmin(rayon1, fmin(rayon2, rayon3));
        double mid = rayonMin + (rayonMax-rayonMin)/2.0;
        printf("RayonMin : %f, RayonMax : %f, Mid : %f \n", rayonMin ,rayonMax, mid);

        // Parcourir tous les noeuds et regarder dans à quel rayon ils appartiennent...

        int* noeudsRayonInterne = malloc(sizeof(int) * theMotorMesh->nElemDomain[Air_gap]*3); // Alloue trop de mémoire mais c'est pas grave
        int* noeudsRayonExterne = malloc(sizeof(int) * theMotorMesh->nElemDomain[Air_gap]*3); // Alloue trop de mémoire mais c'est pas grave
        int indexRayonInterne = 0;
        int indexRayonExterne = 0;


        // marked[i] = 0 si le ieme noeud n'a pas encore été visité et 1 sinon
        int* marked = calloc( theMotorMesh->nNode, sizeof(int)); // Alloue bcp trop de mémoire mais bon
        int count = 0; // nombre de noeuds différents dans Air_Gap
        double rayonNode = 0.0; 

        for(int iTriangle = startTriangleAirGap; iTriangle < endTriangleAirGap; iTriangle++)
        // Parcourir tous les triangles dans Air_Gap 
        {
            for(int i = 0; i < 3; i++)
            // Pour chaque noeud, calculer son rayon
            {
                int iNode = elem[iTriangle*3+i];
                //printf("iNode : %d \n", iNode);
                if( !marked[iNode])
                // le noeud n'a pas déjà été visité
                {
                    //printf("iNode : %d \n", iNode);
                    marked[iNode] = 1;
                    rayonNode = radius(iNode); 
                    //printf("current radius : %f | bound : %f \n",rayonNode,  mid);
                    count++;

                    if( rayonNode > mid) 
                    // noeud situé sur le rayon extérieur du Air_gap
                    {
                        noeudsRayonExterne[indexRayonExterne] = iNode;
                        // ajouter le noeud dans la structure contenant les noeuds sur le rayon externe
                        indexRayonExterne++;
                    }
                    else if(rayonNode < mid )
                    // noeud situé sur le rayon intérieur du Air_gap
                    {
                        noeudsRayonInterne[indexRayonInterne] = iNode;
                        // ajouter le noeud dans la structure contenant les noeuds sur le rayon interne
                        indexRayonInterne++;
                    }
                }
            }
        }
        free(marked);

        printf("Noeuds sur le rayon interne : %d \n", indexRayonInterne);
        printf("Noeuds sur le rayon externe : %d \n", indexRayonExterne);
        printf("Nombre total de noeud : %d \n", count);

        noeudsRayonInterne = realloc((void*) noeudsRayonInterne, sizeof(int) * indexRayonInterne);
        noeudsRayonExterne = realloc((void*) noeudsRayonExterne, sizeof(int) * indexRayonExterne);

        
        printIntArray("Noeuds sur le rayon interne", noeudsRayonInterne);
        for(int i=0; i< indexRayonInterne; i++)
        {
            printf("%d - noeud %d : (%f,%f)\n",i, noeudsRayonInterne[i], X[noeudsRayonInterne[i]], Y[noeudsRayonInterne[i]]);
        }

        printIntArray("Noeuds sur le rayon externe", noeudsRayonExterne);
        for(int i=0; i< indexRayonExterne; i++)
        {
            printf("%d - noeud %d : (%f,%f)\n",i, noeudsRayonExterne[i], X[noeudsRayonExterne[i]], Y[noeudsRayonExterne[i]]);
        }

        // Pour chaque noeud dans le rayon intérieur, déterminer les 2 noeuds les plus proches 
        // dans le rayon extérieur
        int closestNodes[2];  // numéros des noeuds les plus proches
        double closestDistances[2] = {-1.0,-1.0};
        double currentDistance;

        int iNewTriangle = startTriangleAirGap;

        for(int internNode = 0; internNode < indexRayonInterne; internNode++)
        {
            closestNodes[0] = -1;   closestDistances[0] = DBL_MAX;
            closestNodes[1] = -1;   closestDistances[1] = DBL_MAX;

            for(int externNode = 0; externNode < indexRayonExterne; externNode ++)
            {
                // calculer les distances entre internNode et externNode
                currentDistance = squarredDistanceBetweenNodes(internNode, externNode); // la distance au carré est utilisée

                // comparer avec les précédentes distances calculées et mettre à jour si besoin
                if(currentDistance < closestDistances[0])
                {
                    closestDistances[1] = closestDistances[0];
                    closestNodes[1] = closestNodes[0];
                    closestNodes[0] = externNode;
                    closestDistances[0] = currentDistance;
                }
                else if(currentDistance < closestDistances[1])
                {
                    closestNodes[1] = externNode;
                    closestDistances[1] = currentDistance;
                }

            }
            // Créer le triangle maintenant que les 2 noeuds les plus proches du noeud internNode sont connus
            elem[iNewTriangle*3] = internNode;
            elem[iNewTriangle*3+1] = closestNodes[0];
            elem[iNewTriangle*3+2] = closestNodes[1];
            iNewTriangle++;
        }
   
        for(int externNode = 0; externNode < indexRayonExterne; externNode++)
        {
            closestNodes[0] = -1;       closestDistances[0] = DBL_MAX;
            closestNodes[1] = -1;       closestDistances[1] = DBL_MAX;

            for(int internNode = 0; internNode < indexRayonInterne; internNode++)
            {
                currentDistance = squarredDistanceBetweenNodes(externNode, internNode);
                if(currentDistance < closestDistances[0])
                {
                    closestNodes[1] = closestNodes[0];
                    closestDistances[1] = closestDistances[0];
                    closestNodes[0] = internNode;
                    closestDistances[0] = currentDistance;
                }
                else if(currentDistance < closestDistances[1])
                {
                    closestNodes[1] = internNode;
                    closestDistances[1] = currentDistance;
                }
            }
            // Créer le triangle maintenant que les 2 noeuds les plus proches du noeud externNode sont connus
            elem[iNewTriangle*3] = externNode;
            elem[iNewTriangle*3+1] = closestNodes[0];
            elem[iNewTriangle*3+2] = closestNodes[1];
            iNewTriangle++;
        }
        
        if(iNewTriangle - startTriangleAirGap != theMotorMesh->nElemDomain[Air_gap])
        {
            printf("ERROR: Nombre de nouveaux triangles créés: %d, (expected: %d)  \n", iNewTriangle - startTriangleAirGap, theMotorMesh->nElemDomain[Air_gap]);
            // je devrais avoir créé autant de triangle que ceux précédemment présent dans le maillage Air_Gap
        }

    free(noeudsRayonInterne);
    free(noeudsRayonExterne);
    
    
    
}

double motorComputeCouple(motor *theMotor)
{
    return 1e-8;
}

void motorComputeCurrent(motor *theMotor)
{
    return;   
}

/** motorComputeMagneticPotential
 * @in 
 * - motor* theMotor = la structure représentant le moteur étudié
 * @out 
 * -> rempli motor->mesh->a par le potentiel magnétique obtenu en résolvant l'équation de Poisson
 *  en tout point du maillage 
 * 
 * NB : inspirée de la solution du devoir 4 disponible sur https://www.youtube.com/watch?v=580gEIVVKe8.
 */
void motorComputeMagneticPotential(motor* theMotor)
{
    motorMesh* theMotorMesh = theMotor->mesh;
    double* js = theMotor->js;
    double* mu = theMotor->mu;
    int* domain = theMotorMesh->domain;
    int nTriangles = theMotorMesh->nElem;
    int nNode = theMotorMesh->nNode;

    femMesh* theFemMesh = motorMeshToFemMeshConverter(theMotorMesh);
    femEdges* theEdges = femEdgesCreate(theFemMesh);
    femFullSystem* theSystem = femFullSystemCreate(nNode);
    femIntegration* theRule = femIntegrationCreate(3,FEM_TRIANGLE); // règle d'intégration de Hammer à 3 points
    femDiscrete* theSpace = femDiscreteCreate(3,FEM_TRIANGLE); // élément triangulaire bilinéaire

 
    double x[3],y[3],phi[3],dphidxsi[3],dphideta[3],dphidx[3],dphidy[3];
    int iTriangle,iInteg,iEdge,i,j,map[3];
    
    
    for (iTriangle = 0; iTriangle < nTriangles; iTriangle++) 
    {
        myFemMeshLocal(theFemMesh, iTriangle, map, x, y); 

        for (iInteg = 0; iInteg < theRule->n; iInteg++) 
        {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) 
            {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }
            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            for (i = 0; i < theSpace->n; i++) 
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) 
            { 
                for(j = 0; j < theSpace->n; j++) 
                {
                    theSystem->A[map[i]][map[j]] += (1.0 / mu[domain[iTriangle]]) * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight; 
                    // La permabilité dépend des sous-domaines. Il faut donc d'abord déterminer à quel sous-domain appartient le ième triangle (iTriangle).
                    // Il faut ensuite diviser par sa perméabilité pour correspondre avec l'équation à résoudre 
                }
            }                                                                                            
            for (i = 0; i < theSpace->n; i++) 
            {
                theSystem->B[map[i]] += js[domain[iTriangle]] * phi[i] * jac * weight; 
            }
        }
    } 

     for (iEdge= 0; iEdge < theEdges->nEdge; iEdge++) 
     {      
        if (theEdges->edges[iEdge].elem[1] < 0) 
        {  
            for (i = 0; i < 2; i++) 
            {
            	int iNode = theEdges->edges[iEdge].node[i];
            	femFullSystemConstrain(theSystem,iNode,0.0);  
            }
        }
    }

    femFullSystemEliminate(theSystem);
    theMotor->a = theSystem->B;

    freeFemMeshConverted(theFemMesh);
} 

//
// ========= Projet à réaliser ===================
//

// Pistes d'améliorations:
// Enlever du number dans la fonction motorMeshToFemMeshConverter
// Pour ajuster les noeuds dans moving_nodes, comme ce sont les domaines 8,9 et 10 on
// peut éviter de passer sur l'intégralité des noeuds mais uniquement sur ceux la.
