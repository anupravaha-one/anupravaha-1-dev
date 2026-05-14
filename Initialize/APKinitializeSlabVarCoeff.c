#include "APKinitializeSlabVarCoeff.h"

/**
This function creates memory in slabs for variable coefficients like density,viscosity, etc and initializes
it to the user entered initialization value.
*/
double *initializeSlabVarCoeff(int p,Zone *z,double constant,int Count,int dir,int *Point_Range)
{
	double *variable;
	int totalcells,k;
	
	if(p==0)
	{
		totalcells = 1;

		variable=(double *)malloc(totalcells*sizeof(double));

		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=constant;
		}
	}
	else
	{
		totalcells = 0;
		if(dir==1||dir==2)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
		if(dir==3||dir==4)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
		if(dir==5||dir==6)
			totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

		variable=(double *)malloc(totalcells*sizeof(double));

		for(k=0;k<totalcells;k++) /*Here constant =1.0, but these slabs get filled with exact values during RefreshVarCoeff.*/
		{
			*(variable+k)=constant;
		}
	}
	return(variable);
}
