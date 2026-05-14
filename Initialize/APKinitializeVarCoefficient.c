#include "APKinitializeVarCoefficient.h"
#define Cell_Dim 3

/**
This function creates necessary memory for coefficients for Variables like density, viscosity, thermal conductivity etc
and initializes the memory to user entered initialization value.
*/

double *initializeVarCoefficient(int p,Zone *z,double constant,int count)
{
	double *variable;
	int totalcells,k;
	
	if(p==0) /*If loop for checking whether memory is to be allocated or not*/
	{
		totalcells=1;
		
		variable=(double*)malloc(totalcells*sizeof(double)); 
		 
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=constant;
		}         
	}
	else /*p=1 case ie variable case*/
	{
		totalcells=1;
		
		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+k+Cell_Dim));
	
		variable=(double*)malloc(totalcells*sizeof(double));
		
		if(variable == NULL)
		{
			printf("Memory allocation not sufficient,Exiting from InitializeVarCoefficient\n");
			exit(1);
		}
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=constant;
		}
	}
	//printf("I am leaving initvar coeff\n");
	return(variable);
}

