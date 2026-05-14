#include "APKinitializeRadiationSlab.h"

/**
This routine which is called in the initializeRadiation() function performs the task of creating and alloting the sufficient
memory for that radiation variable in Zone 'slab' structure for both Boundary and Interface.
*/

double *initializeRadiationSlab(int i, Base *CgnsBase, int dir, int *Point_Range)
{
	double *variable;
	int totalcells, k, nt, np;
	nt=CgnsBase->BaseParameters->ntheta;
	np=CgnsBase->BaseParameters->nphi;
	if(i==1)
	{
		if(dir==1||dir==2)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1))*nt*np;
		if(dir==3||dir==4)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0))*nt*np;
		if(dir==5||dir==6)
			totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0))*nt*np;
		variable=(double *)malloc(totalcells*sizeof(double));
		if(variable==NULL)
		{  
			printf("\nMemory Insufficient, memory notcreated, exiting from APKinitializeslab\n");
			exit(1);
		}
		
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=0.5;
		}
	}
	else
		variable=NULL;
	
	return(variable);
}
