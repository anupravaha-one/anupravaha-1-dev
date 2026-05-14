#include "APKinitializeRadiationVarCoefficient.h"

/**
This function creates necessary memory for refractive index and initializes the memory to user entered initialization value.
*/


double *initializeRadiationVarCoefficient(int flagvar, Base *CgnsBase, Zone *z, double refindex)
	{
	int totalcells, k, Cell_Dim;
	double *variable;
	//Cell_Dim=CgnsBase->Cell_Dimension;
	//nt=CgnsBase->BaseParameters->ntheta;
	//np=CgnsBase->BaseParameters->nphi;
	//printf("Variable=%e\n",refindex);
	Cell_Dim=3;
		if(flagvar==0)
		{
		totalcells=1;
		variable=(double *)malloc(totalcells*sizeof(double));
			for(k=0;k<totalcells;k++)
			*(variable+k)=refindex;
		}
		else
		{
		 totalcells=1;
		 for(k=0;k<Cell_Dim;k++)
		 	totalcells=totalcells*(*(z->Size+Cell_Dim+k));
		variable=(double *)malloc(totalcells*sizeof(double));
		for(k=0;k<totalcells;k++)
			*(variable+k)=1.0;
		}
	return(variable);
	 }
