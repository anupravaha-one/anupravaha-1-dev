#include "APKinitializeRadiation_variable.h"

/**
This function creates memory for Intensity and initializes them to user entered initialization value.
*/

double *initializeRadiation_variable(int i, Base *CgnsBase, Zone *z, double intensity)
	{
	 double *variable;
	 int k, totalcells, nt, np, Cell_Dim;
	 Cell_Dim=3;//Cell_Dim=CgnsBase->Cell_Dimension;
	 nt=CgnsBase->BaseParameters->ntheta;
	 np=CgnsBase->BaseParameters->nphi;
	 
	 if(i==1) 
		{
		totalcells=1*nt*np;
		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+Cell_Dim+k));
	
		variable=(double*)malloc(totalcells*sizeof(double));
		
		if(variable==NULL)
		{
			printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
			exit(1);
		}
		 
		for(k=0;k<totalcells;k++)
		*(variable+k)=intensity;
		         
		}
	else
		variable=NULL;
		
	return(variable);
	}
