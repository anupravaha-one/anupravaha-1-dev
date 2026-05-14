#include "APKinitializeSlabDispersionCoeff.h"

/**
This function creates memory in slabs for dispersiom coeffcients etc and initializes
it to the user entered initialization value.
*/
double *initializeSlabDispersionCoeff(int p,Zone *z,double *a,int Count,int dir,int *Point_Range)
{
	double *variable;
	int totalcells,k,c;
	
	
	{
		totalcells = 0;
		if(dir==1||dir==2)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+4)- *(Point_Range+1));
		if(dir==3||dir==4)
			totalcells=(*(Point_Range+5)- *(Point_Range+2)) * (*(Point_Range+3)- *(Point_Range+0));
		if(dir==5||dir==6)
			totalcells=(*(Point_Range+4)- *(Point_Range+1)) * (*(Point_Range+3)- *(Point_Range+0));

			
		variable=(double *)malloc(9*totalcells*sizeof(double));

		for(k=0;k<totalcells;k++)
		{       
		         for(c=0;c<9;c++)
			*(variable+9*k+c)=*(a+c);
			
			
		}
	}
	return(variable);
}





double *initializeslab_Delta(int i,int c,int dir,int *Point_Range)
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
		
		
		totalcells= totalcells*c;
		variable=(double *)malloc(totalcells*sizeof(double));
		
		if(variable==NULL)
		{  
			printf("\nMemory Insufficient, memory notcreated, exiting from APKinitializeslab\n");
			exit(1);
		}
		
		
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=0.0;
		
		}
		
		
	}
	else
		variable=NULL;
	
	return(variable);
}
