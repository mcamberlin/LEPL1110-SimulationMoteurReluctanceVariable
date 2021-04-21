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

} 

//
// ========= Projet à réaliser ===================
//