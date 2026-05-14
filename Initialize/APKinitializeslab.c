#include "APKinitializeslab.h"

/**
This routine which is called in the initialize function performs the task of creating and alloting the sufficient
memory for that particular variable in Zone 'slab' structure for both Boundary and Interface.
*/
double *initializeslab(int i,int dir,int *Point_Range)
{
	double *variable;
	int totalcells,k;
	if(i==1)
	{
		
		if(dir==1||dir==2)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
		if(dir==3||dir==4)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
		if(dir==5||dir==6)
			totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

            	if(totalcells < 0)
		  totalcells = -totalcells;

		
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

double *initSlab(int i,int dir,int *Point_Range,double constant)
{
	double *variable;
	int totalcells,k;
	if(i==1)
	{
		
		if(dir==1||dir==2)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
		if(dir==3||dir==4)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
		if(dir==5||dir==6)
			totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

		if(totalcells < 0)
		  totalcells = -totalcells;

		
		variable=(double *)malloc(totalcells*sizeof(double));
		//	printf("\ndirection = %d and total cells = %d\n",dir,totalcells);

		//	printf("\nprange  = [%d %d %d] -> [%d %d %d]\n",*(Point_Range),*(Point_Range+1),*(Point_Range+2),*(Point_Range+3),*(Point_Range+4),*(Point_Range+5));
		
		if(variable==NULL)
		{
			printf("\nMemory Insufficient, memory notcreated, exiting from APKinitializeslab\n");
			exit(1);
		}
		
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=0*constant;
		}
	}
	else
		variable=NULL;
	
	return(variable);
}
