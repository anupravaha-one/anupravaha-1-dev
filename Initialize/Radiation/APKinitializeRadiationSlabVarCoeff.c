#include "APKinitializeRadiationSlabVarCoeff.h"
/**
This function creates memory in slabs for refractive index and initializes it to the user entered initialization value.
*/

double *initializeRadiationSlabVarCoeff(int flagvar, Base *CgnsBase, Zone *z, double refindex, int dir, int *point_Range)
	{
	 int totalcells, k;
	 double *variable;
	 //nt=CgnsBase->BaseParameters->ntheta;
	 //np=CgnsBase->BaseParameters->nphi;
	 if(flagvar==0)
	 	{
		totalcells=1;
		variable=(double *)malloc(totalcells*sizeof(double));
			for(k=0;k<totalcells;k++)
				*(variable+k)=refindex;
		}
	else
		{
		if(dir==1||dir==2)
			totalcells=(*(point_Range+5)-*(point_Range+2))-(*(point_Range+4)-*(point_Range+1));
		if(dir==3||dir==4)	
			totalcells=(*(point_Range+5)-*(point_Range+2))-(*(point_Range+3)-*(point_Range+0));
		if(dir==5||dir==6)	
			totalcells=(*(point_Range+4)-*(point_Range+1))-(*(point_Range+3)-*(point_Range+0));
			
			variable=(double *)malloc(totalcells*sizeof(double));
			
			for(k=0;k<totalcells;k++)
			*(variable+k)=refindex;
		}
	return(variable);
	}
