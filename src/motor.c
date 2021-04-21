#include "motor.h"

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
}

void printFemMesh(femMesh* theFemMesh)
{
    int max = 10;
    printf("\n ====== femMesh informations ============================\n");
    printf("    elem = [ ");
    for(int i =0 ;i < max; i++)
    {
        printf("%d, ",theFemMesh->elem[i]); // PROBLEME ICI
    }
    printf("... ]\n");


    printf("    X = [ ");
    for(int i =0 ;i < max; i++)
    {
        printf("%f, ",theFemMesh->X[i]);
    }
    printf("... ]\n");

    printf("    Y = [ ");
    for(int i =0 ;i < max; i++)
    {
        printf("%f, ",theFemMesh->Y[i]);
    }
    printf("... ]\n");
    printf("    Number of elements           : %d\n",theFemMesh->nElem);   
    printf("    Number of nodes              : %d\n",theFemMesh->nNode);  
    printf("    Number of local nodes        : %d\n",theFemMesh->nLocalNode); 
        
    printf("    number = [ ");
    for(int i =0 ;i < max; i++)
    {
        printf("%d, ",theFemMesh->number[i]);
    }
    printf("... ]\n");
    printf("=========================================================================\n");

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
    double* a = theMotor->a;
    double* js = theMotor->js;
    double* mu = theMotor->mu;
    int nTriangles = theMotorMesh->nElem;

    femMesh* theFemMesh = motorMeshToFemMeshConverter(theMotorMesh);
    
    printFemMesh(theFemMesh);

    femEdges* theEdges = femEdgesCreate(theFemMesh);

    printf("debug\n");

    femFullSystem* theSystem = femFullSystemCreate(theFemMesh->nNode);
    femIntegration* theRule = femIntegrationCreate(3,FEM_TRIANGLE); // règle d'intégration de Hammer à 3 points
    femDiscrete* theSpace = femDiscreteCreate(3,FEM_TRIANGLE); // élément triangulaire bilinéaire

 
    double x[3],y[3],phi[3],dphidxsi[3],dphideta[3],dphidx[3],dphidy[3];
    int iTriangle,iInteg,iEdge,i,j,map[3];
    
    
    for (iTriangle = 0; iTriangle < nTriangles; iTriangle++) 
    {
        femMeshLocal(theFemMesh, iTriangle, map, x, y); 

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
                    theSystem->A[map[i]][map[j]] += (dphidx[i] * dphidx[j] 
                                                   + dphidy[i] * dphidy[j]) * jac * weight; // CHANGEMENT NECESSAIRES ?
                }
            }                                                                                            
            for (i = 0; i < theSpace->n; i++) 
            {
                theSystem->B[map[i]] += phi[i] * jac *weight; // changement nécessaires
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
} 

//
// ========= Projet à réaliser ===================
//