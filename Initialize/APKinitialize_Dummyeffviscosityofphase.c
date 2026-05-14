#include "APKinitialize_Dummyeffviscosityofphase.h"
#define Cell_Dim 3

/**
This function creates memory for Dummy effective viscosity of phase(product of effective viscosity of the phase and 
volume fraction of the phase) which is used as diffusion coefficient for phasic momentum equation in case of turbulent MPFV
and also caluculats the same.
*/

double *initialize_Dummyeffviscosityofphase(int p,Zone *z,double constant,int nophase)
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
		//printf("The value of volume fraction is %lf\n\n",*(z->Flowsol->VarthetaPhase[nophase]));
		flagviscoeff = 0;//for constant viscosity of the phases
		for(k=0;k<NK-1;k++)
		for(j=0;j<NJ-1;j++)
		for(i=0;i<NI-1;i++)
		{
			cellindex = (NI-1)*(NJ-1)*k+(NI-1)*j+i;
			*(variable+cellindex) =*(z->Flowsol->VarthetaPhase[nophase]+cellindex)*(constant + *(z->Flowsol->ViscosityofPhase[nophase]+flagviscoeff*cellindex));
		}
		//printf("the variable is %lf",*(variable));
	}
	else
	{

		variable=NULL;
	}

	return(variable);
}
