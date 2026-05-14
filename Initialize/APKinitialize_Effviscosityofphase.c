#include "APKinitialize_Effviscosityofphase.h"
#define Cell_Dim 3

/**
This function creates memory for Effective viscosity of phase and also calculates the value of mu_q + mu_tq and stores 
in effective viscosity of phase in case of Turbulent MPFV.
 */

double *initialize_Effviscosityofphase(int p,Zone *z,double constant,int nophase)
{
	double *variable;
	int totalcells,NI,NJ,NK,i,j,k,cellindex,flagviscoeff;
	NI = *(z->Size);
	NJ = *(z->Size+1);
	NK = *(z->Size+2);
	if(p==1) /*If loop for checking whether memory is to be allocated or not*/
	{
		totalcells=1;

		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+k+Cell_Dim));

		variable=(double*)malloc(totalcells*sizeof(double));

		if(variable==NULL)
		{
			printf("Memory not Sufficient and exiting from initialize_Effviscosity\n");
			exit(1);
		}
		
		flagviscoeff = 0;//for constant viscosity of the phases
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(variable+cellindex) = constant + *(z->Flowsol->ViscosityofPhase[nophase]+flagviscoeff*cellindex);
		}
	}
	else
	{

		variable=NULL;
	}

	return(variable);
}
