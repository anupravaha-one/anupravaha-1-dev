#include "APKinitializeDispersionCoeff.h"
#define Cell_Dim 3

/**
This function creates necessary memory for coefficients for dispersion coeffceints for anisotropic 
and initializes the memory to user entered initialization value.
*/

double *initializeDispersionCoeff(int p,Zone *z,double *a,int count)
{
	double *variable;
	int totalcells,k,c;
	
		totalcells=1;
		
		for(k=0;k<Cell_Dim;k++)
		totalcells=totalcells*(*(z->Size+k+Cell_Dim));
	
		variable=(double*)malloc(9*totalcells*sizeof(double));
		
		if(variable == NULL)
		{
			printf("Memory allocation not sufficient,Exiting from InitializeVarCoefficient\n");
			exit(1);
		}
		
		for(k=0;k<totalcells;k++)
		{       
		         for(c=0;c<9;c++)
			*(variable+9*k+c)=*(a+c);
		}
	
	//printf("I am leaving initvar coeff\n");
	return(variable);
}

