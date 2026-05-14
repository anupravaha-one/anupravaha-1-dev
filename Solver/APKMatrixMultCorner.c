#include "APKMatrixMultCorner.h"

/**
This function calculates,
LHS = (coefficient_phiP)x(variable_phiP)
RHS = sum of 18 ((coefficient_phiNeighbours)x(variable_phiNeighbours))
which can be used as,
Ax = LHS + RHS
which is the left hand side of any equation Ax = B
This is called only for corner, edge and slab cells.
*/

double MatrixMultCorner(int index[3],int *size,double u_p,double u[18],double *dia,double *DiaN[18],double *LHS,double *RHS)
{
//	printf("entered DiffusionReduced\n");
	int NK,NJ,NI,k,j,i,indexP;
//	double APART;
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	k=index[0],j=index[1],i=index[2];

	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
//	APART = ((*(DiaN[6]+indexP))*u[6]+(*(DiaN[7]+indexP))*u[7]+(*(DiaN[8]+indexP))*u[8]+(*(DiaN[9]+indexP))*u[9]+(*(DiaN[10]+indexP))*u[10]+(*(DiaN[11]+indexP))*u[11]+(*(DiaN[12]+indexP))*u[12]+(*(DiaN[13]+indexP))*u[13]+(*(DiaN[14]+indexP))*u[14]+(*(DiaN[15]+indexP))*u[15]+(*(DiaN[16]+indexP))*u[16]+(*(DiaN[17]+indexP))*u[17]);
	
	*LHS=(*(dia+indexP))*u_p;
	
	*RHS=(*(DiaN[0]+indexP))*u[0]+(*(DiaN[1]+indexP))*u[1]+(*(DiaN[2]+indexP))*u[2]+(*(DiaN[3]+indexP))*u[3]+(*(DiaN[4]+indexP))*u[4]+(*(DiaN[5]+indexP))*u[5];//+APART;
	
	return (0);
}

