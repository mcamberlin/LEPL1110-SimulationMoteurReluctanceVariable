#include "motor.h"

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

static motorMesh* theGlobalMotorMesh;

//
// ========= Fonctions utiles ===================
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

int startAirGap(const motorMesh* theMotorMesh)
{
    int start = 0;
    for(int i=0; i< Air_gap; i++)
    {
        start += theMotorMesh->nElemDomain[i];
    }
    return start; 
}

int endAirGap(const motorMesh* theMotorMesh)
{
    return startAirGap(theMotorMesh) + theMotorMesh->nElemDomain[Air_gap];    
}

/** haveCommonEdge
 * @in 
 * - const motorMesh* theMotorMesh = ensemble des maillages du moteur considéré
 * - const int iTriangle1 = numéro du 1er triangle à considérer
 * - const int iTriangle1 = numéro du 2eme triangle à considérer
 * @out 
 * - retourne un tableau des 2 noeuds formant l'arête commune si elle existe
 * - retourne NULL sinon
 */
int* haveCommonEdge(const motorMesh* theMotorMesh, const int iTriangle1, const int iTriangle2)
{
    int* rslt = malloc( sizeof(int) *2);
    rslt[0] = -1; rslt[1] = -1;

    int* elem = theMotorMesh->elem;
    for(int i1 = 0; i1<3; i1++)
    {
        for(int i2 = 0; i2<3; i2++)
        {
            if(iTriangle1 == 706)
            {
                if(elem[iTriangle1*3 + i1] == elem[iTriangle2*3 + i2])
                {
                    if(rslt[0] == -1)
                    {
                        rslt[0] = elem[iTriangle1*3 + i1];
                    }
                    else
                    {
                        int tmp = rslt[0];
                        rslt[0] = fmin(tmp,elem[iTriangle1*3 + i1]);
                        rslt[1] = fmax(tmp, elem[iTriangle1*3 + i1]);
                        return rslt;
                    }
                }
            }
        }

    }
    return rslt; // aucune edge commune trouvée
}

double distanceBetweenNodes(const motorMesh* theMotorMesh, const int node1, const int node2)
{
    double* X = theMotorMesh->X;
    double* Y = theMotorMesh->Y;
    int* elem = theMotorMesh->elem;

    double x1 = X[elem[node1]];
    double y1 = Y[elem[node1]];

    double x2 = X[elem[node2]];
    double y2 = Y[elem[node2]];

    return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double radius(const motorMesh* theMotorMesh, const int node)
{
    double x = theMotorMesh->X[node];
    double y = theMotorMesh->Y[node];
    return sqrt(x*x + y*y);
}

double angle(const int node)
{
    double* X  = theGlobalMotorMesh->X;
    double* Y = theGlobalMotorMesh->Y;
    int* elem = theGlobalMotorMesh->elem;

    return atan2(Y[elem[node]], X[elem[node]]);
}

void printFemMesh(const femMesh* theFemMesh)
{
    int start = 0;
    int stop = 10;
    printf("\n ====== femMesh informations ============================\n");
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


/**
 * selon l'angle
*/
int compareTo (const void * this, const void * that)
{
    int* thisInt = (int*) this;
    int* thatInt = (int*) that;

    int thisNode = *(thisInt);
    int thatNode = *(thatInt);

    double angleThis = angle(thisNode);
    double angleThat = angle(thatNode);

    return angleThis - angleThat;
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

    // 1. Déterminer les noeuds sur le rayon intérieur ou sur le rayon extérieur
        int startTriangleAirGap = startAirGap(theMotorMesh);
        int endTriangleAirGap = endAirGap(theMotorMesh);
        //printf("startAirGap : %d, endAirGap : %d \n", startTriangleAirGap, endTriangleAirGap);

        printTriangle(theMotorMesh,startTriangleAirGap);

        // calculer de la largeur annulaire : calculer le rayon de 3 noeuds d'un triangle dans Air_gap
        
        double rayon1 = radius(theMotorMesh, elem[startTriangleAirGap*3]);
        double rayon2 = radius(theMotorMesh, elem[startTriangleAirGap*3+1]);
        double rayon3 = radius(theMotorMesh, elem[startTriangleAirGap*3+2]);
        double rayonMax = fmax(rayon1, fmax(rayon2, rayon3));
        double rayonMin = fmin(rayon1, fmin(rayon2,rayon3));
        double mid = rayonMin + (rayonMax-rayonMin)/2.0;
        //printf("RayonMin : %f, RayonMax : %f, Mid : %f \n", rayonMin ,rayonMax, mid);

        int* noeudsRayonInterne = malloc(sizeof(int) * theMotorMesh->nElemDomain[Air_gap]*3); // Alloue trop de mémoire mais c'est pas grave
        int iRayonInterne = 0;
        int* noeudsRayonExterne = malloc(sizeof(int) * theMotorMesh->nElemDomain[Air_gap]*3); // Alloue trop de mémoire mais c'est pas grave
        int iRayonExterne = 0;

        double rayonNode = 0.0;
        for(int iTriangle = startTriangleAirGap; iTriangle < endTriangleAirGap; iTriangle++)
        // Parcourir tous les noeuds dans Air_Gap et calculer leur rayon
        {
            for(int i=0; i< 3; i++)
            {
                rayonNode = radius(theMotorMesh, elem[iTriangle*3+i]);
                //printf("current radius : %f | bound : %f \n",rayonNode,  mid);
                if(rayonNode < mid)
                // noeud situé sur le rayon intérieur du Air_gap
                {
                    noeudsRayonInterne[iRayonInterne] = elem[iTriangle*3]; // ajouter le noeud dans la structure contenant les noeuds sur le rayon interne
                    iRayonInterne++;
                }
                else
                // noeud situé sur le rayon extérieur du Air_gap
                {
                    noeudsRayonExterne[iRayonExterne] = elem[iTriangle*3]; // ajouter le noeud dans la structure contenant les noeuds sur le rayon externe
                    iRayonExterne++;
                }
            }
        }
        noeudsRayonInterne = realloc(noeudsRayonInterne, sizeof(int) * iRayonInterne);
        noeudsRayonExterne = realloc(noeudsRayonExterne, sizeof(int) * iRayonExterne);
        printf("nombre de noeuds sur le rayon interne et externe : %d %d \n", iRayonInterne, iRayonExterne);
        
    // 2. Trier les 2 tableaux noeudsRayonInterne et noeudsRayonExterne selon leur angle
        qsort( (void*) noeudsRayonInterne, iRayonInterne, sizeof(int), &compareTo);
        qsort( (void*) noeudsRayonExterne, iRayonExterne, sizeof(int), &compareTo);


    // 3.
        // Pour chaque noeud dans le rayon exterieur, déterminer les 2 noeuds les plus proches 
        // dans le rayon intérieur
        for(int iNode = 0; iNode < iRayonExterne; iNode++)
        {
            // recherche dichotomique  via l'angle du noeud courant dans noeudsRayonInterne pour trouver  
            // le noeud le plus proche
            int lo = 0;
            int hi = iRayonInterne;
            while(lo <= hi)
            {
                int mid = lo + (hi-lo)/2.0;
                if(angle(iNode) < angle(noeudsRayonInterne[mid])) {hi = mid-1;}
                else if(angle(iNode) > angle(noeudsRayonInterne[mid])) {lo = mid+1;}
                else {}
            }
            // keep track of last compared node to find the closest (when the equals does not exist)
        }

        // Pour chaque noeud dans le rayon intérieur, déterminer les 2 noeuds les plus proches 
        // dans le rayon extérieur
        for(int iNode = 0; iNode < iRayonExterne; iNode++)
        {
            // recherche dichotomique selon l'angle du noeud courant dans noeudsRayonExterne pour trouver
            // le noeud le plus proche
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
// Trier les noeuds dans air_gap pour avoir plus facile pour les retrouver par la suite.


/**
 *  Idée pour trouver l'arête commune:
 * Passer sur tous les noeuds de air_gap 
 * creer une structure triangle a remailler contenant:
 * - les 3 noeuds dans l'ordre croissant
 * - un hash qui correspondant à la combinaison de ces 3 noeuds
 * et trier ses structures selon leur hash.
 */