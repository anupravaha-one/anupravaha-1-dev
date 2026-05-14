#include "APKinitialize_variable.h"
#define Cell_Dim 3

/**
This function creates memory for varibles(u,v,w,p,T,etc) and initializes them to user entered initialization value.
*/
double *initialize_variable(int i,Zone *z,double constant)
{
	double *variable;
	int k,totalcells=1;
	
	if(i==1) /*If loop for checking whether memory is to be allocated or not*/
	{
		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+k+Cell_Dim));
	
		variable=(double*)malloc(totalcells*sizeof(double));
		
		if(variable==NULL)
		{
			printf("Memory Allocation Not Sufficient,Exiting from initialize_variable\n");
			exit(1);
		}
		 
		for(k=0;k<totalcells;k++)
		{
			*(variable+k)=constant;
		}         
	}
	else
	{
		
		variable=NULL;
	}
	
	return(variable);
}
/**
This function creates memory for mass flux.
*/
double *initializeFlux(int i,Zone *z)
{
	double *variable;
	int k,totalcells;
	
	if(i==1) /*If loop for checking whether memory is to be allocated or not*/
	{
		totalcells=1;
		
		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+k+Cell_Dim));
	
		variable=(double*)malloc(6*totalcells*sizeof(double)); 
		
		if(variable==NULL)
		{
			printf("Memory Allocation Not Sufficient,Exiting from InitializeFlux in initialize_variable.c\n");
			exit(1);
		}
		 
		for(k=0;k<(6*totalcells);k++)
		{
			*(variable+k)=0.0;
		}
	}
	else
	{
		variable=NULL;
	}
	
	return(variable);
}

//Added by sushil for diffusion flux
double *initializeDiffFlux(int i,Zone *z)
{
	double *variable;
	int k,totalcells;
// 	nop=z->Nphase;
	if(i==1) /*If loop for checking whether memory is to be allocated or not*/
	{
		totalcells=1;
		
		for(k=0;k<Cell_Dim;k++)
			totalcells=totalcells*(*(z->Size+k+Cell_Dim));
	
		variable=(double*)malloc(totalcells*sizeof(double)); 
		
		if(variable==NULL)
		{
			printf("Memory Allocation Not Sufficient,Exiting from InitializeDiffFlux in initialize_variable.c\n");
			exit(1);
		}
		 
		for(k=0;k<(totalcells);k++)
		{
			*(variable+k)=0.0;
		}
	}
	else
	{
		variable=NULL;
	}
	
	return(variable);
}


