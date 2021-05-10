#include "motor.h"
#include <float.h>

typedef int bool;
#define true 1
#define false 0

#define Stator_core 0
#define Rotor_core 8   
#define Coil_AP 1
#define Coil_AN 2
#define Coil_BP 3
#define Coil_BN 4
#define Coil_CP 5
#define Coil_CN 6
#define Rotor_gap 10
#define Air_gap 11
#define Stator_gap 7
#define Rotor_air 9


const static motorMesh* theGlobalMotorMesh;

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
        theFemMesh->number = malloc(sizeof(int) * theFemMesh->nNode); 
        for (int i = 0; i < theFemMesh->nNode; i++)
        {
            theFemMesh->number[i] = i; 
        } 
        return theFemMesh;
    }

    void freeFemMeshConverted(femMesh* theFemMesh)
    {
        //free(theFemMesh->number);
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

    // NOTE : 
    // Par triangle intérieur, on désigne un triangle dont 2 de ses noeuds sont situés sur le rayon intérieur de la bande glissante
    // Par triangle extérieur, on désigne un triangle dont 2 de ses noeuds sont situés sur le rayon extérieur de la bande glissante

    /* triangle
    = structure pour décrire les triangles intérieurs constructibles autour du noeud intérieur node
    */
    typedef struct {
        int node;            // numéro du noeud intérieur
        int* indexTriangles; // indices des triangles auquel appartient le noeud intérieur node
        int nbElems;         // nombre d'éléments dans indexTriangles
    } triangle;

    static triangle** theGlobalNewTriangles;

    /** createNewTriangle
     * @in 
     * - int iNode = numéro du noeud intérieur à ajouter dans la structure
     * - int iTriangle = numéro du triangle auquel le iNode appartient car il est le plus proche
     * @ out
     * - triangle* = pointeur de structure triangle dans lequel le noeud iNode a été ajouté, de la mémoire a été allouée pour 3 autres index dans
     *  indexTriangles, le nombre d'éléments a été initialisé à 1.
     */
    triangle* createNewTriangle(int iNode, int iTriangle)
    {
        triangle* theTriangle = malloc(sizeof(triangle));
        theTriangle->node = iNode;
        theTriangle->indexTriangles = malloc(sizeof(int)*3); // un noeud appartient au plus à 2 autres triangles intérieur (mais j'ai alloué 3 au cas où)
        theTriangle->indexTriangles[0] = iTriangle;
        theTriangle->indexTriangles[1] = -1;
        theTriangle->indexTriangles[2] = -1;
        theTriangle->nbElems = 1;
        return theTriangle;
    }  

    /** freeNewTriangle
     * @in 
     * - triangle* t = pointeur vers la structure triangle à libérer
     * @ out
     * -> désalloue la mémoire occupée par t
     */
    void freeNewTriangle(triangle* t)
    {
        free(t->indexTriangles);
        free(t);
    }

    /** startAirGap
     * @in 
     * @out
     * -retourne l'indice du tableau elem à partir duquel le domaine Air_gap commence
     */
    int startAirGap()
    {
        int start = 0;
        for(int i=0; i< Air_gap; i++)
        {
            start += theGlobalMotorMesh->nElemDomain[i];
        }
        return start; 
    }

    /** endAirGap
     * @in 
     * @out
     * -retourne l'indice +1 du tableau elem à partir duquel le domaine Air_gap termine 
     */
    int endAirGap()
    {
        return startAirGap() + theGlobalMotorMesh->nElemDomain[Air_gap];    
    }

    /** radius
     * @in
     * - const int node = numéro d'un noeud étudié
     * @out
     * - retourne le rayon du noeud
     */
    double radius(const int node)
    {
        double x = theGlobalMotorMesh->X[node];
        double y = theGlobalMotorMesh->Y[node];
        return sqrt(x*x + y*y);
    }

    // valeur de décision pour déterminer si le noeud est sur le rayon intérieur ou sur le rayon extérieur
    static double theGlobalThreshold; 

    /** computeThreshold
     * @in 
     * - const int triangleAirGap = numéro d'un triangle situé dans le domaine Air_gap
     * @out
     * -> met dans la variable statique @theGlobalThreshold le seuil à partir duquel on détermine si un
     * noeud est situé:
     * soit entre le domaine Air_gap et le domaine Stator_gap 
     * soit entre le domaine Air_gap et le domaine Rotor_gap 
     */
    void computeThreshold(const int triangleAirGap)
    {
        int* elem = theGlobalMotorMesh->elem;
        double rayon1 = radius(elem[triangleAirGap*3]);
        double rayon2 = radius(elem[triangleAirGap*3+1]);
        double rayon3 = radius(elem[triangleAirGap*3+2]);
        double rayonMax = fmax(rayon1, fmax(rayon2, rayon3));
        double rayonMin = fmin(rayon1, fmin(rayon2, rayon3));
        theGlobalThreshold = rayonMin + (rayonMax-rayonMin)/2.0;
        //printf("RayonMin : %f, RayonMax : %f, Threshold : %f \n", rayonMin ,rayonMax, theGlobalThreshold);
    }

    /** isExternTriangle
     * @in 
     * - const int iTriangle = numéro du triangle considéré
     * - int* noeudsExt = pointeur vers un tableaux de 2 éléments qui contiendra les numéros des noeuds sur le rayon extérieur
     * - int* indexNoeudsExt = pointeur vers un tableaux de 2 éléments qui contiendra les indexes dans le triangle iTrangle des 2 noeuds sur le rayon extérieur
     * - int* noeudsInt = pointeur vers un tableaux de 2 éléments qui contiendra les numéros des noeuds sur le rayon intérieur
     * - int* indexNoeudsInt = pointeur vers un tableaux de 2 éléments qui contiendra les indexes dans le triangle iTrangle des 2 noeuds sur le rayon intérieur
     * @out 
     * -> met à jour les tableaux passés en arguments
     * et
     * - retourne 1 si le triangle est formé de 2 noeuds sur le rayon extérieur
     * - retourne 0 si le triangle est formé de 2 noeuds sur le rayon intérieur
     */
    int isExternTriangle(const int iTriangle, int* noeudsExt, int* indexNoeudsExt, int* noeudsInt, int* indexNoeudsInt)
    {
            int* elem = theGlobalMotorMesh->elem;

            // Réinitialiser les noeuds courants sur le rayon intérieur/extérieur
            noeudsExt[0] = -1;          noeudsInt[0] = -1;
            noeudsExt[1] = -1;          noeudsInt[1] = -1;
            indexNoeudsExt[0] = -1;     indexNoeudsInt[0] = -1;
            indexNoeudsExt[1] = -1;     indexNoeudsInt[1] = -1;
            double nbNoeudsExt = 0;     double nbNoeudsInt = 0;

            int iNode; 
            double rayonNode;
            
            for(int i=0; i<3; i++)
            {
                iNode = elem[iTriangle*3+i]; 
                //printf("iNode : %d \n", iNode);

                rayonNode = radius(iNode); 
                //printf("current radius : %f | bound : %f \n",rayonNode,  theGlobalThreshold);

                if(rayonNode > theGlobalThreshold) // noeud situé sur le rayon extérieur du Air_gap
                {
                    nbNoeudsExt++;
                    if(noeudsExt[0] == -1)
                    {
                        noeudsExt[0] = iNode;
                        indexNoeudsExt[0] = i;
                    }
                    else if(noeudsExt[1] == -1)
                    {
                        noeudsExt[1] = iNode;
                        indexNoeudsExt[1] = i;
                    }
                }
                else if(rayonNode < theGlobalThreshold) // noeud situé sur le rayon intérieur du Air_gap
                {
                    nbNoeudsInt++;
                    if(noeudsInt[0] == -1)
                    {
                        noeudsInt[0] = iNode;
                        indexNoeudsInt[0] = i;
                    }
                    else if(noeudsInt[1] == -1)
                    {
                        noeudsInt[1] = iNode;
                        indexNoeudsInt[1] = i;
                    }
                }
            }

            if( nbNoeudsExt == 2 && nbNoeudsInt == 1)
            {
                return 1;
            }
            else if( nbNoeudsExt == 1 && nbNoeudsInt == 2)
            {
                return 0;
            }
            else
            // should never happen
            {
                return -1;
            }
    }

    /** distance
     * @in
     * const double x1, y1 = coordonnées du point 1
     * const double x2, y2 = coordonnées du point 2
     * @out
     * - retourne la distance entre les 2 points (x1,y1) et (x2,y2)
     */
    double distance(const double x1, const double y1, const double x2, const double y2)
    {
        return sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    }

    /** squarredDistance
     * @in
     * const double x1, y1 = coordonnées du point 1
     * const double x2, y2 = coordonnées du point 2
     * @out
     * - retourne la distance au carré entre les 2 points (x1,y1) et (x2,y2)
     */
    double squarredDistance(const double x1, const double y1, const double x2, const double y2)
    {
        return  (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
    }

//
// ========= Fonctions utiles motorComputeCouple() ===================
//

    #define MU_0 12.5663706147*1e-7 //[kg m A−2 s−2] 
    // source: https://fr.wikipedia.org/wiki/Constante_magnétique

    /** startRotorGap
     * @in 
     * @out
     * -retourne l'indice du tableau elem à partir duquel le domaine Rotor_gap commence
     */
    int startRotorGap()
    {
        int start = 0;
        for(int i=0; i< Rotor_gap; i++)
        {
            start += theGlobalMotorMesh->nElemDomain[i];
        }
        return start; 
    }

    /** endRotorGap
     * @in 
     * @out
     * -retourne l'indice +1 du tableau elem à partir duquel le domaine Rotor_gap termine 
     */
    int endRotorGap()
    {
        return startRotorGap() + theGlobalMotorMesh->nElemDomain[Rotor_gap];    
    }


    /** computeWidth
     * @in 
     * - const int triangleAirGap = numéro d'un triangle situé dans le domaine Rotor_gap
     * @out
     * - retourne la largeur de la bande Rotor_gap
     */
    double computeWidth(const int triangleRotorGap)
    {
        int* elem = theGlobalMotorMesh->elem;
        double rayon1 = radius(elem[triangleRotorGap*3]);
        double rayon2 = radius(elem[triangleRotorGap*3+1]);
        double rayon3 = radius(elem[triangleRotorGap*3+2]);
        double rayonMax = fmax(rayon1, fmax(rayon2, rayon3));
        double rayonMin = fmin(rayon1, fmin(rayon2, rayon3));
        //printf("RayonMin : %f, RayonMax : %f, Largeur : %f \n", rayonMin ,rayonMax, rayonMax-rayonMin);
        return rayonMax-rayonMin;
    }


//
// ========= Fonctions utiles motorComputeCoupleBandSolver() ===================
//
    femDiffusionProblem * myFemDiffusionCreate(femMesh* theFemMesh)
    {    
        femDiffusionProblem *theProblem = malloc(sizeof(femDiffusionProblem));
        theProblem->mesh  = theFemMesh;
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);        // utilisation de triangle
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE);      // règle de Hammer à 3 points
        theProblem->size = theFemMesh->nNode;
        theProblem->sizeLoc = theProblem->mesh->nLocalNode;
        femMeshRenumber(theProblem->mesh,FEM_XNUM);                 // renumérotation selon X
        theProblem->sourceValue = 1.0;                              
        theProblem->dirichletValue = 0.0;                           // valeur imposée sur la frontière est de 0
        theProblem->dirichlet = calloc(theProblem->size, sizeof(int));    
        femEdges *theEdges = femEdgesCreate(theProblem->mesh);
        for (int i = 0; i < theEdges->nEdge; i++) 
        {      
            if (theEdges->edges[i].elem[1] < 0) 
            {       
                theProblem->dirichlet[theEdges->edges[i].node[0]] = 1; 
                theProblem->dirichlet[theEdges->edges[i].node[1]] = 1; 
            }
        }
        femEdgesFree(theEdges);
        int band = femMeshComputeBand(theProblem->mesh);
        theProblem->solver = femSolverBandCreate(theProblem->size, theProblem->sizeLoc,band); 
        theProblem->soluce = calloc(theProblem->size,sizeof(double));
        return theProblem;
    }


//
// ========= Projet à réaliser ===================
//

void motorAdaptMesh(motor *theMotor, double delta)
{
    motorMesh *theMotorMesh = theMotor->mesh;  
    theGlobalMotorMesh = theMotorMesh;  // mets à jour la variable statique
    
    double* X = theMotorMesh->X;
    double* Y = theMotorMesh->Y;
    int* elem = theMotorMesh->elem;
    int* nElemDomain = theMotorMesh->nElemDomain;
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
    
    int startTriangleAirGap = startAirGap();
    int endTriangleAirGap = endAirGap();
    //printf("startAirGap : %d, endAirGap : %d \n", startTriangleAirGap, endTriangleAirGap);
    //printTriangle(theMotorMesh,startTriangleAirGap);
    
    // calculer de la valeur seuil en calculant le rayon de 3 noeuds d'un triangle dans Air_gap
    computeThreshold(startTriangleAirGap);

    int noeudsExt[2];           int indexNoeudsExt[2];
    int noeudsInt[2];           int indexNoeudsInt[2];


    double xOpt; // abscisse optimale pour créer un triangle équilatéral
    double yOpt; // ordonnée optimale pour créer un triangle équilatéral

    int iNode2;
    double currentDistance;
    double minDistance; // distance la plus proche entre un noeud intérieur et extérieur
    int closestIntNode; // noeud sur le rayon intérieur le plus proche

    triangle** theNewTriangles = malloc(sizeof(triangle*) * nElemDomain[Air_gap]);
    theGlobalNewTriangles = theNewTriangles;
    int nbNewTriangles = 0;

    for(int iTriangle = startTriangleAirGap; iTriangle < endTriangleAirGap; iTriangle++)
    // Parcourir tous les triangles dans Air_Gap à la recherche de triangle formé de
    // 2 noeuds sur le rayon extérieur
    {
    
        if( isExternTriangle(iTriangle, noeudsExt, indexNoeudsExt, noeudsInt, indexNoeudsInt))
        {
            // Déterminer la position optimale du noeud sur le rayon intérieur pour former un triangle équilatéral:
            xOpt = (X[noeudsExt[0]] + X[noeudsExt[1]]) /2.0;
            yOpt = (Y[noeudsExt[0]] + Y[noeudsExt[1]]) /2.0;

            // Trouver le noeud sur le rayon interne qui est le plus proche du point optimal   
            minDistance = DBL_MAX;
            closestIntNode = -1; 
            for(int iTriangle2 = startTriangleAirGap; iTriangle2 < endTriangleAirGap; iTriangle2++)
            {
                for(int i = 0; i<3; i++)
                {
                    iNode2 = elem[iTriangle2*3+i];
                    if(radius(iNode2) < theGlobalThreshold)
                    // Si le noeud se situe sur le rayon intérieur
                    {
                        // calculer les distances entre le point optimal et le noeud courant
                        currentDistance = squarredDistance(xOpt, yOpt, X[iNode2], Y[iNode2]);
                        if(currentDistance < minDistance)
                        // comparer avec les précédentes distances calculées et mettre à jour si besoin
                        {
                            minDistance = currentDistance;
                            closestIntNode = iNode2;
                        }
                    }
                }
            }  
            
            // Nouveau triangle trouvé formé des noeuds: noeudsExt[0], noeudsExt[1] et closestIntNode
            //printf("iTriangle %d: noeuds (%d %d %d) \t après rotation: (%d %d %d) \n",iTriangle, noeudsExt[0], noeudsExt[1], noeudInt, noeudsExt[0], noeudsExt[1], closestIntNode);
            elem[iTriangle*3+indexNoeudsInt[0]] = closestIntNode;

            
            // Déterminer si closestNode est déja ds un triangle
            bool isInNewTriangles = false;
            int indexInNewTriangles = -1;
            for(int i = 0; i < nbNewTriangles; i++)
            {
                if(theNewTriangles[i]->node == closestIntNode)
                {
                    isInNewTriangles = true;
                    indexInNewTriangles = i;
                    break;
                }
            }
            
            if(! isInNewTriangles)
            // closestNode n'est pas dans theNewTriangles
            {
                triangle* newTriangle = createNewTriangle(closestIntNode, iTriangle);
                theNewTriangles[nbNewTriangles] = newTriangle;
                nbNewTriangles++;
            }
            
            else
            // le noeud optimal est déjà dans theNewTriangles
            {
                triangle* tmp = theNewTriangles[indexInNewTriangles];
                int* tmpIndexTriangles = tmp->indexTriangles;
                if(tmp->nbElems > 2)
                {
                    printf("ERROR: la structure triangle ne contient pas assez de place pour tous les triangles constructibles à partir de ce noeud\n");
                }
                tmpIndexTriangles[ tmp->nbElems] = iTriangle;
                tmp->nbElems++;
            }        
        }
    }


    int noeudsExtTriangleIntern[2];        int noeudsIntTriangleIntern[2]; 
    int indexNoeudsExtTriangleIntern[2];   int indexNoeudsIntTriangleIntern[2];

    int indexNodeInt1InTheNewTriangles;
    int indexNodeInt2InTheNewTriangles;
    int nbMatch = 0; // variable utilisée pour éviter de parcourir tous les éléments de theNewTriangle

    for(int iTriangle = startTriangleAirGap; iTriangle < endTriangleAirGap; iTriangle++)
    // Reparcourir tous les triangles dans Air_Gap à la recherche de triangle intérieur
    {
        if( ! isExternTriangle(iTriangle, noeudsExtTriangleIntern, indexNoeudsExtTriangleIntern, noeudsIntTriangleIntern, indexNoeudsIntTriangleIntern))
        // Pour un triangle intérieur
        {
            indexNodeInt1InTheNewTriangles = -1;
            indexNodeInt2InTheNewTriangles = -1;
            nbMatch = 0;

            // Parcourir tous les éléments de theNewTriangle pour trouver les positions ds theNewTriangle 
            // des 2 noeuds sur le rayon intérieur
            for(int i=0; i< nbNewTriangles && nbMatch < 2; i++)
            {
                if(theNewTriangles[i]->node == noeudsIntTriangleIntern[0])
                {
                    indexNodeInt1InTheNewTriangles = i;
                    nbMatch++;
                }
                else if(theNewTriangles[i]->node == noeudsIntTriangleIntern[1])
                {
                    indexNodeInt2InTheNewTriangles = i;
                    nbMatch++;
                }
            }

            // Normalement chaque noeud intérieur doit être présent dans exactement 2 triangles
            if(nbMatch != 2)
            {
                printf("ERROR - Chaque noeud intérieur doit être présent dans exactement 2 triangles MAIS présent dans: %d \n", nbMatch);
            }

            triangle* newTriangle1 = theNewTriangles[indexNodeInt1InTheNewTriangles];
            int* indexTrianglesNewTriangle1 = newTriangle1->indexTriangles;
            int nbElemensNewTriangle1 = newTriangle1->nbElems;

            for(int iElemNewTriangle1 = 0; iElemNewTriangle1 < nbElemensNewTriangle1; iElemNewTriangle1++)
            // Parcourir tous les autres triangles auxquels appartient le noeud NodeInt1InTheNewTriangles
            {
                int iTriangleNewTriangle1 = indexTrianglesNewTriangle1[iElemNewTriangle1];               

                // Déterminer quels sont les noeuds sur le rayon extérieur
                for(int i=0; i< 3; i++)
                {
                    int iNodeNewTriangle1 = elem[iTriangleNewTriangle1 *3 + i];
                    if( radius(iNodeNewTriangle1) > theGlobalThreshold)
                    // si le noeud est sur le rayon extérieur
                    {
                        triangle* newTriangle2 = theNewTriangles[indexNodeInt2InTheNewTriangles];
                        int* indexTrianglesNewTriangle2 = newTriangle2->indexTriangles;
                        int nbElemensNewTriangle2 = newTriangle2->nbElems;

                        for(int iElemNewTriangle2 = 0; iElemNewTriangle2 < nbElemensNewTriangle2; iElemNewTriangle2++)
                        // Parcourir tous les autres triangles auxquels appartient le noeud NodeInt2InTheNewTriangles
                        {
                            int iTriangleNewTriangle2 = indexTrianglesNewTriangle2[iElemNewTriangle2];

                            // Déterminer quels sont les noeuds sur le rayon extérieur
                            for(int j =0; j< 3; j++)
                            {
                                int iNodeNewTriangle2 = elem[iTriangleNewTriangle2 *3 + j];
                                if(  iNodeNewTriangle1 ==  iNodeNewTriangle2 && radius(iNodeNewTriangle2) > theGlobalThreshold)
                                {
                                    elem[iTriangle*3 + indexNoeudsExtTriangleIntern[0]] = iNodeNewTriangle1;    
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // libérer la mémoire allouée dans theNewTriangles
    for(int i =0; i<nbNewTriangles;i++)
    {
        freeNewTriangle(theNewTriangles[i]);
    }
    free(theNewTriangles);   
}

/** motorComputeCouple
 * inspirée de:
 * - la solution du devoir 4 Poisson disponible sur https://www.youtube.com/watch?v=580gEIVVKe8
 * - la solution du devoir 2 Inertial Gear 
 */
double motorComputeCouple(motor *theMotor)
{
    double L = theMotor->L;   
    motorMesh* theMotorMesh = theMotor->mesh;
    theGlobalMotorMesh = theMotorMesh;
    
    double* a = theMotor->a;
    int* elem = theMotorMesh->elem;
        
    int startTriangleRotorGap = startRotorGap();
    int endTriangleRotorGap = endRotorGap();

    double d = computeWidth(startTriangleRotorGap);

    femMesh* theFemMesh = motorMeshToFemMeshConverter(theMotorMesh);
    femIntegration* theRule = femIntegrationCreate(3,FEM_TRIANGLE);   // règle d'intégration de Hammer à 3 points
    femDiscrete* theSpace = femDiscreteCreate(3,FEM_TRIANGLE);          // élément triangulaire bilinéaire

 
    double x[3];
    double y[3];
    int map[3];
    double Jacobien;

    double xsi;    
    double eta;    
    double weight;
    
    double phi[3],dphidxsi[3],dphideta[3],dphidx[3],dphidy[3];
    
    double xLoc[theRule->n]; // Abscises des points d'intégration
    double yLoc[theRule->n]; // Oordonnées des points d'intégration
    double weightsAjusted[theRule->n];
    
    double dxdxsi;          double dydxsi;
    double dxdeta;          double dydeta;
    double jac;

    double r;           double theta;
    double dxdr;        double dxdtheta;  
    double dydr;        double dydtheta;
    double dadr;        double dadtheta;

    double couple = 0.0;

    for (int iTriangle = startTriangleRotorGap; iTriangle < endTriangleRotorGap; iTriangle++) 
    // Parcourir tous les éléments dans Air_Gap
    {
        myFemMeshLocal(theFemMesh, iTriangle, map, x, y); 
        Jacobien = fabs( (x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]) ); 
        // J = abs( (X2-X1).(Y3-Y1) - (Y2-Y1).(X3-X1) )

        for (int iInteg = 0; iInteg < theRule->n; iInteg++) 
        {    

            xsi    = theRule->xsi[iInteg];
            eta    = theRule->eta[iInteg];
            weight = theRule->weight[iInteg]; 

            //Points d'intégrations
            xLoc[iInteg] = x[0] * (1- xsi - eta) + x[1] * xsi + x[2] * eta;
            yLoc[iInteg] = y[0] * (1- xsi - eta) + y[1] * xsi + y[2] * eta;
            //Poids ajustés
            weightsAjusted[iInteg] = weight * Jacobien;
            //Calcul intermédiaire de l'intégrale
            //I += weightsAjusted[iInteg] * f(xLoc[iInteg], yLoc[iInteg]);
 
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            dxdxsi = 0.0;         dxdeta = 0.0;
            dydxsi = 0.0;         dydeta = 0.0;
            for (int i = 0; i < theSpace->n; i++) 
            {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
            }

            jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (int i = 0; i < theSpace->n; i++) 
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }   
            
            r = sqrt(xLoc[iInteg]*xLoc[iInteg] + yLoc[iInteg]*yLoc[iInteg]);
                                        theta = atan2(yLoc[iInteg],xLoc[iInteg]);
            dadr = 0.0;                 dadtheta = 0.0;
            for (int i = 0; i < theSpace->n; i++) 
            { 
                dadtheta += a[map[i]] * (dphidx[i] * -r*sin(theta) + dphidy[i] * r*cos(theta));
                dadr += a[map[i]] * (dphidx[i] * cos(theta) + dphidy[i] * sin(theta));
            }           
            couple +=  dadr * dadtheta * weightsAjusted[iInteg];                           
        }
    } 
    couple = (-L/(MU_0 * d) ) *couple ; 
    //printf("Couple : %f [Nm]\n", couple);
    return couple;
}


void motorComputeCurrent(motor *theMotor)
{
    double js = 8.8464*1e5;                         // densité de courant [A/m**2]
    double angle = fmod( theMotor->theta, 2*M_PI); // angle électrique du rotor [rad]
    //printf("angle : %f° \n", angle*180/M_PI);
    angle = fmod( theMotor->theta, M_PI);
    angle = fabs(angle);

    if (0.0 <= angle && angle < M_PI/6.0)
    // [0°, 30°[
    {
        theMotor->js[Coil_AP] = -js;
        theMotor->js[Coil_AN] = js;
        theMotor->js[Coil_BP] = js;
        theMotor->js[Coil_BN] = -js;
        theMotor->js[Coil_CP] = 0;
        theMotor->js[Coil_CN] = 0;
    } 
    else if (M_PI/6.0 <= angle && angle < 2*M_PI/6.0)
    // [30°, 60°[
    {
        theMotor->js[Coil_AP] = -js;
        theMotor->js[Coil_AN] = js;
        theMotor->js[Coil_BP] = 0.0;
        theMotor->js[Coil_BN] = 0.0;
        theMotor->js[Coil_CP] = js;
        theMotor->js[Coil_CN] = -js;
    } 
    else if (2*M_PI/6.0 <= angle && angle < M_PI/2.0)
    // [60°, 90°[
    {
        theMotor->js[Coil_AP] = 0.0;
        theMotor->js[Coil_AN] = 0.0;
        theMotor->js[Coil_BP] = -js;
        theMotor->js[Coil_BN] = js;
        theMotor->js[Coil_CP] = js;
        theMotor->js[Coil_CN] = -js;
    }
    else if (M_PI/2.0 <= angle && angle < 4*M_PI/6.0)
    // [90°, 120°[
    {
        theMotor->js[Coil_AP] = js;
        theMotor->js[Coil_AN] = -js;
        theMotor->js[Coil_BP] = -js;
        theMotor->js[Coil_BN] = js;
        theMotor->js[Coil_CP] = 0;
        theMotor->js[Coil_CN] = 0;
    } 
    else if (4*M_PI/6.0 <= angle && angle < 5*M_PI/6.0)
    // [120°, 150°[
    {
        theMotor->js[Coil_AP] = js;
        theMotor->js[Coil_AN] = -js;
        theMotor->js[Coil_BP] = 0;
        theMotor->js[Coil_BN] = 0;
        theMotor->js[Coil_CP] = -js;
        theMotor->js[Coil_CN] = js;
    } 
    else if (5*M_PI/6.0 <= angle && angle < M_PI)
    // [150°, 180°[
    {
        theMotor->js[Coil_AP] = 0;
        theMotor->js[Coil_AN] = 0;
        theMotor->js[Coil_BP] = js;
        theMotor->js[Coil_BN] = -js;
        theMotor->js[Coil_CP] = -js;
        theMotor->js[Coil_CN] = js;
    }    
    return;   
}

/** motorComputeMagneticPotentialFullSolver
 * inspirée de la solution du devoir 4 Poisson disponible sur https://www.youtube.com/watch?v=580gEIVVKe8.
 * utilisation d'une élimination gaussienne

void motorComputeMagneticPotentialFullSolver(motor* theMotor)
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
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
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

    // libérer toutes les mémoires allouées:
    freeFemMeshConverted(theFemMesh);
    
    //femEdgesFree(theEdges);
    //femFullSystemFree(theSystem);
    //femIntegrationFree(theRule);
    //femDiscreteFree(theSpace);      
}   

*/
/** motorComputeMagneticPotential
 * inspirée de la solution du devoir 5 BandSolver
 * utilisation d'un solver bande
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
    femDiffusionProblem* theProblem = myFemDiffusionCreate(theFemMesh);
    //femSolverPrintInfos(theProblem->solver);

    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femSolver *theSolver = theProblem->solver;

    int *number = theProblem->mesh->number;
    double source = theProblem->sourceValue;
    double dirichlet = theProblem->dirichletValue;

 
    double phi[3],dphidxsi[3],dphideta[3],dphidx[3],dphidy[3];
    double Xloc[3],Yloc[3], Uloc[3];
    int iTriangle,iInteg,iEdge,i,j,map[3],ctr[3];
    
    
    double **A = theSolver->local->A;
    double *Aloc = theSolver->local->A[0];
    double *Bloc = theSolver->local->B;
    
    for (iTriangle = 0; iTriangle < theMotorMesh->nElem; iTriangle++) 
    {
        for (i = 0; i < theSpace->n; i++)
        {
            Bloc[i] = 0;
        }      
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++)
        {
            Aloc[i] = 0;
        } 
        
        femDiffusionMeshLocal(theProblem,iTriangle,map,ctr,Xloc,Yloc,Uloc);  
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) 
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
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; 
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) 
            {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) 
            { 
                for(j = 0; j < theSpace->n; j++) 
                {
                    A[i][j] += (1.0 / mu[domain[iTriangle]]) * (dphidx[i] * dphidx[j] + dphidy[i] * dphidy[j]) * jac * weight; 
                    // La permabilité dépend des sous-domaines. Il faut donc d'abord déterminer à quel sous-domain appartient le ième triangle (iTriangle).
                    // Il faut ensuite diviser par sa perméabilité pour correspondre avec l'équation à résoudre 
                }
            }                                                                                            
            for (i = 0; i < theSpace->n; i++) 
            {
                Bloc[i] += js[domain[iTriangle]] * phi[i] * jac * source * weight; 
            }
        }
        for (i = 0; i < theSpace->n; i++) 
        {
            if (ctr[i] == 1) 
            {
                femFullSystemConstrain(theSolver->local,i,dirichlet);
            }
        }    
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); 
    } 
 
    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->size; i++)
    {
        theProblem->soluce[i] += soluce[number[i]];
        theMotor->a[i] = theProblem->soluce[i];
    }
    
    // libérer toutes les mémoires allouées:
    freeFemMeshConverted(theFemMesh);
    //femDiffusionFree(theProblem);
}

void motorFree(motor *theMotor)
{
    free(theMotor->mesh->elem);
    free(theMotor->mesh->X);
    free(theMotor->mesh->Y);
    free(theMotor->mesh->nElemDomain);
    free(theMotor->mesh->nameDomain);
    free(theMotor->mesh->domain);
    free(theMotor->mesh);
    free(theMotor->a);
    free(theMotor->movingNodes);
    free(theMotor->js);
    free(theMotor->mu);
    free((void*) theMotor->hystereticCurveH);
    free((void*) theMotor->hystereticCurveB);
    free(theMotor);
    return;
}