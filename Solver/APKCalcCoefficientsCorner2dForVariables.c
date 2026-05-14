#include "APKCalcCoefficientsCorner2dForVariables.h"

/**
This function calculates coeffients of diffusion term in the convection-diffusion equation for velocities, temperature
and scalars. This is used for the case of constant diffusion coefficients. It is used for 2D case. It is called
only for corner, edge and slab cells.
Note:  It just multiplies the diffusion coefficient (e.g., kinematic viscosity) with the pressure Poisson equation
coefficients which are calculated by CallCalcCoefficients().
*/

double CalcCoefficientsCorner2dForVariables(Zone *zone,int index[3],double vol_p,double u[18],int *size,double *beta,double *DiaN[18],int vartag)
{
	int NK,NJ,NI,k,j,i,indexP,l;
	double *diffusionCoeffdomain,*CoeffN[18],*diffdia,diag;
	
	NI=*(size+0);
	NJ=*(size+1);
	NK=*(size+2);
	
	diffusionCoeffdomain = MakeDiffusionCoeffDomain(zone,vartag);
	
	diffdia = (zone->grid)->DiffusionCoeff_P;
	for(l=0;l<18;l++)
	{
		CoeffN[l] = (zone->grid)->DiffusionCoeff_nb[l];
	}
	
	k=index[0],j=index[1],i=index[2];
	
	indexP=(NI-1)*(NJ-1)*k+(NI-1)*j+i;
	
	diag=(*diffusionCoeffdomain)*(*(diffdia+indexP));
		 
	*(DiaN[0]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[0]+indexP));
		     
	*(DiaN[1]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[1]+indexP));
		     
	*(DiaN[2]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[2]+indexP));
		     
	*(DiaN[3]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[3]+indexP));
		      
	*(DiaN[4]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[4]+indexP));
	
	*(DiaN[5]+indexP)=(*diffusionCoeffdomain)*(*(CoeffN[5]+indexP));
	
	*(DiaN[6]+indexP)=0.0;
	
	*(DiaN[7]+indexP)=0.0;
	
	*(DiaN[8]+indexP)=0.0;
	
	*(DiaN[9]+indexP)=0.0;
	
	*(DiaN[10]+indexP)=0.0;
	
	*(DiaN[11]+indexP)=0.0;
	
	*(DiaN[12]+indexP)=0.0;
	
	*(DiaN[13]+indexP)=0.0;
	
	*(DiaN[14]+indexP)=0.0;
	
	*(DiaN[15]+indexP)=0.0;
	
	*(DiaN[16]+indexP)=0.0;
	
	*(DiaN[17]+indexP)=0.0;
	
	return (diag);
}


